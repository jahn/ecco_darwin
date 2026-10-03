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
#   mixing          GGL90 (not KPP) -- read from the archive, see code_offline_gmredi
#   irradiance      pkg/oasim computes spectral surface irradiance at runtime;
#                   pkg/radtrans propagates it through 13 wavebands, 400-700 nm
#
# Derived from Steph Dutkiewicz's 1-degree offline Darwin configuration, ported to
# LLC90 and to real ECCOv4r6 physics over the full 1992-2025 record (her original
# was a repeating 360-day climatology from ECCO Iteration 3.73).
#
# Directory contents:
#
#   code_v4r6_offline      compile-time headers, packages.conf, NAS optfiles
#   code_offline_gmredi    pkg/offline patched to read generic vertical diffusivity
#                          and all GMRedi tensor components and bolus velocity
#   input_v4r6_offline/    namelists + PBS job script for the offline Darwin run
#   input_v4r6_forward/    data.diagnostics for the ECCOv4r6 forward run -- the one
#                          file that differs from the official release (see Part 4)
#   util/regrid/           regrids the 1-degree and LLC270 initial conditions onto
#                          LLC90 (run once, off-line, before Part 7)
#   util/offline           scripts for preparing offline forcing fields from online run
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

# Then replace the official data.diagnostics with this repo's copy -- it is the
# same 168 official streams plus 18 daily streams (indices 169-186) that the
# offline run needs.  Nothing else in code/ or namelist/ changes.

cp <ecco_darwin>/offline/input_v4r6_forward/data.diagnostics namelist/data.diagnostics


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
# Only input_v4r6_forward/data.diagnostics.  It adds eighteen daily-mean streams
# at indices 169-186 -- THETA, SALT, UVELMASS, VVELMASS, WVELMASS, GGL90Kr,
# GM_Kwx/Kwy/Kwz/Kux/Kvy/Kuz/Kvz/PsiX/PsiY, SIarea, oceFWflx -- writing to
# diags/D*.  Together these are the complete forcing set pkg/offline and
# pkg/darwin need.
#
# SIarea is taken from this run rather than from a separate sea-ice product
# because ECCOv4r6 already runs pkg/seaice: the archive gives real, native-LLC90,
# daily ice cover over the whole record.  The official release has SIarea only as
# a monthly stream (index 60), which is not sufficient here.  oceFWflx is the
# freshwater flux into/out of the ocean needed to reproduce the changes in seasurface
# height that in turn drive the dynamic stretching of the model grid.
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
# compact.
# Only four fields need processing: horizontal velocities since the offline
# package expects mass transports divided by the at-rest grid stretching factors
# h0FacW and h0FacS, and seaice area and freshwater flux as these are read in
# using exf and need to be assembled into yearly files.
# Yearly files are preferred over one single large file for all years as this
# cycling over a single, repeated year of physical state easier.

cd util/offline
module load python3/3.9.12
pip install --user MITgcmutils          # once
export V4R6_RUN_DIR=<WORKINGDIR>/MITgcm/ECCOV4/release6

# assemble daily averages into yearly fields for exf
for f in DSIarea DoceFWflx; do ./mkexfdaily.py $V4R6_RUN_DIR/diags/$f {1992..2025}; done

# compute offline velocities from mass transports: uveloff=uvelmass/h0FacW
./mkofflineveldaily.py $V4R6_RUN_DIR 24 298056 24

# link first daily average (iteration 24) to previous day (may want to link last too)
./link24 $V4R6_RUN_DIR/diags


# ========
# Part 6 -- Build the offline Darwin executable
#
# This build uses darwin3, not stock MITgcm.  The two builds are independent --
# they share no MITgcm checkpoint, because pkg/offline communicates through plain
# binary arrays rather than checkpoint-specific pickup formats.
# With an offline run, there generally isn't much reason to run with the same checkpoint
# or backport.  Therea are however hacks to stretching factors and the shelfice package
# that would require a bit of attention when switching to current darwin3.
# We have to run the shelfice package though much of it could be disabled in offline mode.
# This just needs to be coded.  So stick with backport_ckpt68g for now.
# The strategy here is to keep any physical options, modification and parameters that
# may still be relevant in offline mode, and only add things from the biogeochemical setup
# that directly affect the biology, so darwin* radtrans* oasim* and header files, and
# the associated runtime files.  The parameter ALLOW_OLD_VIRTUALFLUX in DARWIN_OPTIONS.h
# depends on the physics and is only used for the _linear_ free surface.

git clone https://github.com/darwinproject/darwin3
cd darwin3/pkg/darwin
git checkout backport_ckpt68y
cd ../../..

mkdir -p build_offline_llc90 run_offline_llc90
cd build_offline_llc90
module load comp-intel/2020.4.304 mpi-hpe/mpt \
            hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
setenv MOD $ECCO_DARWIN_OFFLINE      # <ecco_darwin>/offline
../darwin3/tools/genmake2 -rootdir=../darwin3 \
    -optfile="$MOD/code_darwin_offline/linux_amd64_ifort+mpi_ice_nas.electra_skylake_20251013" \
    -mods="$MOD/code_offline_gmredi $MOD/code_v4r6_offline" -mpi
make depend
make -j 16

# -rootdir is required; genmake2 cannot auto-detect ROOTDIR from the relative
# ../darwin3/tools/genmake2 invocation path.
#
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
cp $MOD/input_v4r6_offline/* .

# these are all from the physics run
ln -sf <ancillary_data_root>/native_grid_files/tile*.mitgrid .
ln -sf <ancillary_data_root>/input_init/* .
ln -sf <ancillary_data_root>/input_forcing/other/*.bin .
ln -sf <ancillary_data_root>/input_forcing/exf/* .
ln -sf <ancillary_data_root>/misc/tools/mkdir_subdir_diags.py .
python mkdir_subdir_diags.py

# we run straight off the diagnostics output of the physics run, plus preprossesed fields there
ln -sf $V4R6_RUN_DIR/diags offline

# 3. initial conditions, regridded onto LLC90 -- see util/regrid/ below.
ln -sf <regridded_llc90>/*.bin .

# 4. iron deposition: Hamilton et al. (2020), on its native 144x96 grid.
#    pkg/darwin interpolates at runtime from data.darwin's iron_lon0/lat0 block,
#    so no regridding is needed.
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/CESM-MIMI_1980-2015_CAM4-6MEAN_MonthlyDep_Hamiltonetal2020_clim.R4_144x96x12.bin .

# 5. same for atmospheric pCO2:
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/apCO2_* .

# 6. and BGC river runoff:
ln -sf /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff .

# 7. and OASIM atmospheric input:
ln -sf /nobackup/ojahn/forcing/oasim .

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
