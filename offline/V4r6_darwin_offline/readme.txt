# Offline ECCO-Darwin: LLC90 biogeochemistry driven by pre-computed ECCOv4r6 physics
#
# Darwin 3 integrated with pkg/offline against a daily-mean physical archive
# written by an ECCOv4r6 forward run, rather than stepping the dynamics itself.
#
#   grid            LLC90, 50 levels (sNx=sNy=30, nPx=113) -- same grid as ECCOv4r6
#   period          1992-01-01 to 2025, gregorian calendar, deltaT = 3600 s,
#                   nenditer = 298031
#   ecosystem       36 passive tracers = 20 biogeochemical + 10 plankton carbon
#                   (c01-c10) + 6 chlorophyll; 6 phytoplankton + 4 zooplankton,
#                   no explicit bacteria
#   mixing          GGL90 (not KPP) -- read from the archive, see code_offline
#   irradiance      pkg/oasim computes spectral surface irradiance at runtime;
#                   pkg/radtrans propagates it through 13 wavebands, 400-700 nm
#
# Derived from Steph Dutkiewicz's 1-degree offline Darwin configuration, ported to
# LLC90 and to real ECCOv4r6 physics over the full 1992-2025 record (her original
# was a repeating 360-day climatology from ECCO Iteration 3.73).
#
# Directory contents:
#
#   code_offline/          pkg/offline patched to read missing GMRedi tensor components,
#                          mass volume source and GGL90 vertical diffusivity
#   code_darwin_offline/   additional code mods for 6p4z setup
#   input_offline/         namelists + PBS job script for the generic offline run
#   input_darwin_offline/  additional namelists for 6p4z Darwin run
#   code_v4r6_forward/     code mods for forward run to write addMass diagnostic
#                          needed for offline shelfice
#   input_v4r6_forward/    data.diagnostics for the ECCOv4r6 forward run
#   util/offline/          convert forward-run diagnostic output for offline run
#                          - convert mass transports to offline velocities
#                          - assemble fresh-water flux and seaice area files into
#                          - yearly files for exf
#   util/regrid/           regrids the 1-degree and LLC270 initial conditions onto
#                          LLC90 (run once, off-line, before Part 7)
#
# Three stages, run in sequence:
#   ECCOv4r6 forward run (Parts 1-4) -> convert (Part 5) -> offline Darwin (Parts 6-7)


# ========
# Part 1 -- Get the code (ECCOv4r6 forward run)

mkdir WORKINGDIR
cd WORKINGDIR

git clone https://github.com/MITgcm/MITgcm.git -b checkpoint68g
git clone https://github.com/MITgcm-contrib/ecco_darwin.git

cd MITgcm
mkdir -p ECCOV4/release6
cd ECCOV4/release6
git clone https://github.com/ECCO-GROUP/ECCO-v4-Configurations.git
mv "ECCO-v4-Configurations/ECCOv4 Release 6/code" .
mv "ECCO-v4-Configurations/ECCOv4 Release 6/namelist" .
rm -rf ECCO-v4-Configurations

# Add code modifications for writing the addMass diagnostics needed to reproduce the
# meltwater flux from shelfice

cp <ecco_darwin>/offline/V4r6_darwin_offline/code_v4r6_forward/* code/

# Then replace the official data.diagnostics with this repo's copy -- it is the
# same 168 official streams plus 11 daily streams (indices 169-179) that the
# offline run needs.  Nothing else in code/ or namelist/ changes.

cp <ecco_darwin>/offline/V4r6_darwin_offline/input_v4r6_forward/data.diagnostics namelist/data.diagnostics


# ========
# Part 2 -- Compile the ECCOv4r6 forward run

module purge
module load comp-intel/2020.4.304
module load mpi-hpe/mpt
module load hdf4/4.2.12
module load hdf5/1.8.18_mpt
module load netcdf/4.4.1.1_mpt
module list

mkdir build
cd build
../../../tools/genmake2 -mods=../code \
    -optfile=../code/linux_amd64_ifort+mpi_ice_nas.electra_skylake_20251013 -mpi
make depend
make all
cd ..

# The Skylake optfile matches "model=sky_ele" in Part 3 and will fail with an
# illegal-instruction error on non-Skylake nodes; use the unsuffixed
# linux_amd64_ifort+mpi_ice_nas instead if Skylake cannot be guaranteed.


# ========
# Part 3 -- Run the ECCOv4r6 forward run
#
# Save as run_script.csh (csh, not bash) in MITgcm/ECCOV4/release6/ and qsub it.
# Runtime is roughly 36-39 h on 113 ranks for the full 1992-2025 record.

#PBS -S /bin/csh
#PBS -l select=3:ncpus=40:model=sky_ele
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -o ./
#PBS -m bea

limit stacksize unlimited
module purge
module load comp-intel/2020.4.304
module load mpi-hpe/mpt
module load hdf4/4.2.12
module load hdf5/1.8.18_mpt
module load netcdf/4.4.1.1_mpt
module load python3/3.9.12
module list

setenv FORT_BUFFERED 1
setenv MPI_BUFS_PER_PROC 128
setenv MPI_DISPLAY_SETTINGS
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HOME}/lib
unsetenv MPI_IB_RECV_MSGS
unsetenv MPI_UD_RECV_MSGS

set nprocs   = 113
set basedir  = ./
set inputdir = /nobackup/owang/runs/V4r6/PO.DAAC/ancillary_data/ancillary_data_orig/
set rundir   = ${basedir}/run

if ( -d ${rundir} ) then
 echo 'Directory "run" exists.'
 echo 'Please rename/remove it and re-submit the job.'
 exit 1
endif

mkdir ${rundir}
cd ${rundir}

ln -s ../namelist/* .
ln -s ${inputdir}/input_init/* .
ln -s ${inputdir}/misc/tools/mkdir_subdir_diags.py .
ln -s ${inputdir}/data_constraints/data_error/*/* .
ln -s ${inputdir}/data_constraints/*/* .
ln -s ${inputdir}/input_forcing/other/*.bin .
ln -s ${inputdir}/input_forcing/control_weights/* .
ln -s ${inputdir}/input_forcing/control_weights/atm_ctrls/* .
ln -s ${inputdir}/input_forcing/exf/* .
ln -s ${inputdir}/native_grid_files/tile*.mitgrid .

python mkdir_subdir_diags.py
mkdir -p diags/archive_{theta,salt,uvel,vvel,wvel,gmkwx,gmkwy,gmkwz,ggl90kr,sflux,siarea}
cp -p ../build/mitgcmuv .
mpiexec -np ${nprocs} /u/scicon/tools/bin/mbind.x ./mitgcmuv

# The explicit "mkdir -p diags/archive_*" is needed because mkdir_subdir_diags.py
# only knows about the official streams.
#
# Leave the data.ctrl / data.exf swaps documented in the official reproduction
# instructions commented out: the defaults apply the existing optimal correction
# in a single forward pass, which is what the offline archive requires.
#
# Verify before proceeding: PBS log ends in NORMAL END, and each of the eleven
# diags/archive_*/ directories holds 12418 daily records, evenly spaced 24
# iterations apart, from iteration 24 to 298032.


# ========
# Part 4 -- What differs from the official ECCOv4r6 recipe
#
# Only input_v4r6_forward/data.diagnostics.  It adds eleven daily-mean streams
# at indices 169-179 -- THETA, SALT, UVEL, VVEL, WVEL, GM_Kwx/Kwy/Kwz, GGL90Kr,
# SFLUX and SIarea -- writing to diags/archive_*.  Together these
# are the complete forcing set pkg/offline and pkg/darwin need.
#
# SIarea is taken from this run rather than from a separate sea-ice product
# because ECCOv4r6 already runs pkg/seaice: the archive gives real, native-LLC90,
# daily ice cover over the whole record.  The official release has SIarea only as
# a monthly stream (index 60), which is not sufficient here.
#
# util/archive/data.diagnostics.snippet holds just these streams, for merging by
# hand if the official data.diagnostics has moved on.
#
# The build, compile, run procedure and every input file are otherwise unmodified
# official ECCOv4r6.


# ========
# Part 5 -- Convert the archive into pkg/offline format
#
# No regridding: the forward run's diagnostic output is already native LLC90
# compact.  The files are already single precision and organized as the offline
# package expects them.  With an rstar coordinate, however, mass transports need
# to be converted to "offline velocities" (U-transport/h0FacW, etc, where h0FacW
# is the stretching factor for a seasurface at rest).  And some files need to be
# assembled for the exf package to be able to read the.
#
# in the forward run directory, do this:

module load python3/3.9.12
pip install --user MITgcmutils          # once

export ECCO_DARWIN_OFFLINE=<ecco_darwin>/offline/V4r6_darwin_offline
export UTIL=$ECCO_DARWIN_OFFLINE/util/offline
$UTIL/mkofflineveldaily.py . 24 298056 24

for f in DSIarea DoceFWflx; do $UTIL/mkexfdaily.py diags/$f {1992..2025}; done

$UTIL/mklink24 diags


# ========
# Part 6 -- Build the offline Darwin executable
#
# This build uses darwin3, not stock MITgcm.  The two builds are independent --
# they share no MITgcm checkpoint, because pkg/offline communicates through plain
# binary arrays rather than checkpoint-specific pickup formats.

git clone https://github.com/darwinproject/darwin3
cd darwin3/pkg/darwin
git checkout backport_ckpt68y
cd ../../..

mkdir -p build_offline_llc90 run_offline_llc90
cd build_offline_llc90
module load comp-intel/2020.4.304 mpi-hpe/mpt \
            hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
setenv MOD $ECCO_DARWIN_OFFLINE      # <ecco_darwin>/offline/V4r6_darwin_offline
../darwin3/tools/genmake2 -rootdir=../darwin3 \
    -optfile="$MOD/code_darwin_offline/linux_amd64_ifort+mpi_ice_nas.electra_skylake_20251013" \
    -mods="$MOD/code_offline $MOD/code_darwin_offline" -mpi
make depend
make -j 16

# On Pleiades, darwin3's own code generator (tools/darwin/cogapp/cogapp.py) calls
# hashlib.md5() at three sites without usedforsecurity=False and dies under
# FIPS-mode OpenSSL ("EVP_DigestInit_ex ... disabled for FIPS").  MD5 is used
# there only as an internal cache key.  Patch before building:
#
#   sed -i 's/hashlib.md5()/hashlib.md5(usedforsecurity=False)/' \
#       darwin3/tools/darwin/cogapp/cogapp.py
#
# This must be reapplied after any fresh darwin3 clone.


# ========
# Part 7 -- Run the offline Darwin simulation

cd ../run_offline_llc90
ln -sf ../build_offline_llc90/mitgcmuv .
cp -p $MOD/input_offline/* .
cp -fp $MOD/input_darwin_offline/* .

# 1. grid and bathymetry.  tile*.mitgrid comes from the same ancillary-data root
#    as Part 3.  data's bathyFile is the ECCOv4r5/r6 LLC90 bathymetry,
#    BATHY_ICE_SHELF_CAVITY_PLUS_ICE_FRONT_LLC_0090.bin -- verified correct for
#    r6 (its RF cumulative depths match r6's own delR exactly).
ln -sf <ancillary_data_root>/native_grid_files/tile*.mitgrid .
ln -sf <llc90_grid_dir>/BATHY_ICE_SHELF_CAVITY_PLUS_ICE_FRONT_LLC_0090.bin .

# 2. the offline fields are in the run directory of the forward run
ln -sfT <rundir>/diags offline

# 3. initial conditions, regridded onto LLC90 -- see util/regrid/ below.
ln -sf <regridded_llc90>/*.bin .

# 4. iron deposition: Hamilton et al. (2020), on its native 144x96 grid.
#    pkg/darwin interpolates at runtime from data.darwin's iron_lon0/lat0 block,
#    so no regridding is needed.
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/CESM-MIMI_1980-2015_CAM4-6MEAN_MonthlyDep_Hamiltonetal2020_clim.R4_144x96x12.bin .

# 5. atmospheric pCO2: NOAA Marine Boundary Layer, native 256x2 grid, also
#    runtime-interpolated.  pCO2File is read through EXF_INTERP_READ, a
#    direct-access reader that needs ONE file holding all records -- the source
#    is one file per calendar year, so concatenate in chronological order.
#    Expect exactly 36653056 bytes (17897 days, 1979-2027, x 2048 bytes/record).
#    data.darwin's pCO2startdate1 is already set to 1979-01-01, this file's own
#    first record.
cat `ls /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/apCO2_[0-9]* | sort -t_ -k2 -n` > apCO2_combined

# 6. BGC river runoff, 10 fields (DOC/DON/DOP/DIN/DIP/DSi/POC/POP/PON/DIC), from
#    ecco_darwin's own v06/1deg configuration -- already LLC90-native despite the
#    directory name.  These have no interpolation block, so darwin_exf_load.F
#    reads them via MDS_READ_FIELD, which needs a <prefix>.data/.meta pair rather
#    than the raw per-year source files.  Build them:
setenv OFFLINE_RUN_DIR `pwd`
python3 $MOD/util/archive/build_runoff_mds.py

# 7. OASIM atmospheric input.  Same problem again: the source is one file per
#    calendar year (and one file per waveband per year for the three spectral
#    aerosol fields), and the shared source directory is read-only, so this
#    builds a real local oasim/ rather than a symlink.
python3 $MOD/util/archive/build_oasim_mds.py

# 8. optical lookup tables (darwin_waterabsorbFile / darwin_phytoabsorbFile /
#    darwin_particleabsorbFile) -- three small wavelength-indexed text tables,
#    13 wavebands 400-700 nm, matching RT_wbRefWLs.  Copy into the run directory.
#    optics_water.txt, optics_plankton.txt, optics_detritus.txt
#
# 9. hydrothermal vent iron (ventHe3File): mitgcm_3he_flux.bin, 128x64 x 12
#    months, float32 (393216 bytes).  Requires DARWIN_ALLOW_HYDROTHERMAL_VENTS,
#    already set in code_darwin_offline/DARWIN_OPTIONS.h.

qsub job_offline_llc90

# Unlike the forward run, this directory is itself the run directory; the job
# script does not create a run/ subdirectory.
#
# walltime is set to 96:00:00.  A 1-year test measured 0.76 s/step on 113 ranks,
# projecting to about 63 h for the full 298031 steps; the remainder is margin for
# checkpoint and diagnostic writes.  Re-profile and tighten once a full run has
# completed.
#
# For a short test, reduce nenditer in "data" (the full-record value is 298031)
# and cut walltime to match.


# ========
# util/regrid -- building the LLC90 initial conditions
#
# Run once, before Part 7, on a machine with the source files and with
# ecco_v4_py + MITgcmutils available.  Output lands in regridded_llc90/, named to
# match the bare filenames in data.ptracers and data.darwin.
#
#   setenv LLC90_GRID_DIR  <ECCOv4 LLC90 grid: XC/YC/RF/hFacC>
#   setenv LLC270_GRID_DIR <LLC270 grid, for the two LLC270-native tracers>
#   python3 run_regrid_steph_ics.py       # 19 tracers + wind + biomass, 1-deg -> LLC90
#   python3 run_regrid_v06_llc270_ics.py  # DIC and ALK only, LLC270 -> LLC90
#
# Two different methods, because the two sources are on different grids:
# bilinear from the 1-degree lat-lon fields, and nearest-wet-neighbour (KDTree
# over Cartesian-embedded lon/lat) for the LLC270-native DIC and ALK.
#
# The remaining tracers restart from an earlier 1-degree Darwin run and from
# GLODAP / Levitus climatology; data.ptracers documents the per-tracer source in
# a trailing comment on each PTRACERS_initialFile line.  Reordering the tracer
# list breaks that mapping silently.
#
# util/regrid/plot_ic_surface_fields.py, plot_ic_profiles.py and
# plot_forcing_timeseries.py check the regridded output before it is used.


# ========
# External inputs, summarised
#
# Nothing below is in this repository.  All paths are on Pleiades.
#
#   ECCOv4r6 ancillary data (grid, bathymetry, forcing, ICs, pickups)
#       /nobackup/owang/runs/V4r6/PO.DAAC/ancillary_data/ancillary_data_orig/
#   iron deposition, Hamilton et al. (2020)
#       /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/
#   OASIM atmospheric fields
#       /nobackup/ojahn/forcing/oasim/{C61extrap/R4,bcs}
#   atmospheric pCO2, NOAA MBL
#       /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/
#   BGC river runoff
#       ecco_darwin v06/1deg configuration
#   1-degree tracer ICs and the three optics tables
#       from Steph Dutkiewicz; regridded ICs are produced by util/regrid/
#   hydrothermal vent 3He flux
#       mitgcm_3he_flux.bin, 128x64x12, float32


# ========
# Coupled parameters -- change these together
#
# Nothing validates these at run time in a useful way; a mismatch appears as a
# crash or, worse, as silently wrong indexing.
#
#   nplank=10       DARWIN_SIZE.h    = sum(grp_nplank) in data.darwin
#                                    = the c01-c10 names in data.ptracers
#   nPhoto=6        DARWIN_SIZE.h    = groups with grp_photo=1 = the Chl01-Chl06 names
#   nGroup=7        DARWIN_SIZE.h    = length of grp_names and every grp_*(:) array
#   PTRACERS_num=36 PTRACERS_SIZE.h  = PTRACERS_numInUse in data.ptracers
#                                    = 20 + nplank + nPhoto
#   nlam=13         RADTRANS_SIZE.h  = count of RT_wbRefWLs and of each
#                                      darwin_scatSlope* list; RT_wbEdges is nlam+1
#   Nr=50           SIZE.h           = length of delR in "data"; RT_kmax must be <= Nr
#   nPx*nPy=113     SIZE.h           = the -np given to mpiexec, in both job scripts


# ========
# Status
#
# Completed and verified: the ECCOv4r6 forward run (Parts 1-3), the archive
# conversion (Part 5, all 11 fields at 12418 records), and the offline build
# (Part 6).  Part 7 has been exercised as a short test and as a full 1-year test,
# both ending normally.  The full 1992-2025 offline integration has not yet been
# run to completion.
