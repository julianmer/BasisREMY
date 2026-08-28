#!/usr/bin/env Rscript
# spant_worker.R — BasisREMY's bridge into spant (R).
#
#   Rscript spant_worker.R job.json out.json
#
# The same script serves both runtimes: the user's own R (Rscript on the host)
# and the Docker image basisremy-spant (rocker/r-ver + spant). Python writes
# a JSON job, this script runs spant::sim_basis for every requested molecule
# and writes the time-domain signals back as JSON — no R<->Python binding to
# compile on the user's machine.
#
# job.json
#   sequence      "press" | "press_shaped" | "steam" | "slaser" | "spin_echo" |
#                 "mega_press" | "pulse_acquire"
#   ft_hz         transmitter frequency [Hz]
#   fs_hz         sampling frequency (spectral width) [Hz]
#   n             number of complex points
#   ref_ppm       ppm reference of the rotating frame (4.65)
#   linewidth_hz  Lorentzian linewidth applied to every spin group (optional)
#   te_s, te1_s, te2_s, te3_s, tm_s      sequence timings [s] (per sequence)
#   steam_variant "ideal" | "cof" | "young"
#   edit_on_ppm, edit_off_ppm, edit_bw_hz, edit_steps   MEGA-PRESS
#   pulse_file, pulse_format ("pta" | "bruker" | "ascii"), pulse_dur_s,
#   refoc_flip_deg                                        PRESS shaped
#   metabolites   spant molecule names (get_mol_names())
#
# out.json
#   {"ok": true, "spant_version": "4.4.0",
#    "basis": {"<name>": {"re": [...], "im": [...], "spant_name": "NAA"}}}
#   MEGA-PRESS: {"<name>": {"on": {"re", "im"}, "off": {"re", "im"}, ...}}
#   {"ok": false, "error": "..."} on failure.

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript spant_worker.R job.json out.json")
job_path <- args[1]
out_path <- args[2]

write_out <- function(x) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null"),
             out_path)
}

run_job <- function(job) {
  suppressPackageStartupMessages(library(spant))

  acq <- def_acq_paras(ft = job$ft_hz, fs = job$fs_hz, N = as.integer(job$n),
                       ref = job$ref_ppm %||% 4.65)

  # ---- molecules (spant's own spin-system library) ------------------------
  requested <- as.character(job$metabolites)
  lw <- job$linewidth_hz
  mols <- lapply(requested, function(name) {
    # lipid / macromolecule definitions take the field into account (ft);
    # the others have no such argument
    mol <- tryCatch(get_mol_paras(name, ft = job$ft_hz),
                    error = function(e) get_mol_paras(name))
    if (!is.null(lw)) {
      for (i in seq_along(mol$spin_groups)) mol$spin_groups[[i]]$lw <- lw
    }
    mol
  })

  # ---- pulse sequence -----------------------------------------------------
  seq_key <- job$sequence
  steam_variant <- job$steam_variant %||% "ideal"
  pul_seq <- switch(seq_key,
    press         = seq_press_ideal,
    press_shaped  = seq_press_2d_shaped,
    steam         = switch(steam_variant,
                           cof   = seq_steam_ideal_cof,
                           young = seq_steam_ideal_young,
                           seq_steam_ideal),
    slaser        = seq_slaser_ideal,
    spin_echo     = seq_spin_echo_ideal,
    mega_press    = seq_mega_press_ideal,
    pulse_acquire = seq_pulse_acquire,
    stop("unknown sequence '", seq_key, "'"))
  seq_args <- switch(seq_key,
    press         = list(TE1 = job$te1_s, TE2 = job$te2_s),
    press_shaped  = list(TE1 = job$te1_s, TE2 = job$te2_s,
                         pulse_file = job$pulse_file,
                         pulse_dur = job$pulse_dur_s,
                         pulse_file_format = job$pulse_format,
                         refoc_flip_angle = job$refoc_flip_deg %||% 180),
    steam         = list(TE = job$te_s, TM = job$tm_s),
    slaser        = list(TE1 = job$te1_s, TE2 = job$te2_s, TE3 = job$te3_s),
    spin_echo     = list(TE = job$te_s),
    mega_press    = list(TE1 = job$te1_s, TE2 = job$te2_s,
                         BW = job$edit_bw_hz %||% 110,
                         steps = as.integer(job$edit_steps %||% 50)),
    pulse_acquire = list())

  # ---- simulate and pull the time-domain signals out --------------------
  simulate <- function(extra) {
    basis <- do.call(sim_basis, c(list(mol_list = mols, pul_seq = pul_seq,
                                       acq_paras = acq), seq_args, extra))
    mrs <- basis2mrs_data(basis)
    if (is_fd(mrs)) mrs <- fd2td(mrs)
    arr <- mrs$data
    d <- dim(arr)
    nd <- length(d)
    # mrs_data arrays are [x, y, z, coil, dynamic, receiver, N]: one basis
    # signal per dynamic — locate that axis by its length to stay robust
    dyn_axis <- which(d[-nd] == length(requested))[1]
    if (is.na(dyn_axis)) dyn_axis <- 5
    fids <- list()
    for (i in seq_along(requested)) {
      idx <- as.list(rep(1, nd))
      idx[[dyn_axis]] <- i
      idx[[nd]] <- seq_len(d[nd])
      fid <- as.vector(do.call(`[`, c(list(arr), idx)))
      fids[[requested[i]]] <- list(re = Re(fid), im = Im(fid),
                                   spant_name = basis$names[i])
    }
    fids
  }

  if (seq_key == "mega_press") {
    on  <- simulate(list(ed_freq = job$edit_on_ppm))
    off <- simulate(list(ed_freq = job$edit_off_ppm))
    basis_out <- setNames(lapply(requested, function(nm) {
      list(on = on[[nm]][c("re", "im")], off = off[[nm]][c("re", "im")],
           spant_name = on[[nm]]$spant_name)
    }), requested)
  } else {
    basis_out <- simulate(list())
  }

  list(ok = TRUE, spant_version = as.character(packageVersion("spant")),
       basis = basis_out)
}

res <- tryCatch({
  job <- jsonlite::fromJSON(job_path, simplifyVector = TRUE)
  run_job(job)
}, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
write_out(res)
