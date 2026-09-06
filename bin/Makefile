#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
#NVECTOR = 32
#NDIM = 3
#NPRE = 8
#NVAR = 16
#NENER = 0
#SOLVER = hydro
#PATCH = ../patch/Horizon5-master
#EXEC = ramses_ch
USE_FFTW ?= 1
CIC_ATOMIC_SCATTER ?= 0
##################################
# Yohan's patch
# [DE] Use coarse vector tiles; NVECTOR=32 made init_tree launch
# [DE] millions of tiny OpenMP regions for the large zoom IC load.
NVECTOR = 500
NDIM = 3
NPRE = 8
# rho + 3 momentum + total energy + total metal + delayed-cooling reservoir
# + 11 tracked elements.  Unused passive fields remain inert in legacy mode.
NVAR = 18
DUST_LIVE ?= 0
ifeq ($(DUST_LIVE),1)
  # Keep the default legacy field count unchanged.  The live profile reserves
  # the complete channel window plus dust mass and dust thermal energy.
  NVAR = 30
endif
NENER = 0
SOLVER = hydro
PATCH = ../patch/lagRamses
EXEC = ramses_final
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git rev-parse HEAD)$(shell git diff-index --quiet HEAD -- || printf '%s' '-dirty')
GITREPO = $(shell git config --get remote.origin.url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) -DNVAR=$(NVAR) \
          -DSOLVER=$(SOLVER)
#############################################################################
# Fortran compiler: always Intel ifx
#   Default:      64-bit Morton key (level <= 22)
#   MORTON128=1:  128-bit Morton key via dual-int8 mkey_t (level <= 43)
#############################################################################
F90 = mpiifx -qopenmp
ifdef FDMDEBUG
  FFLAGS = -fpp -I../patch/lagRamses -O2 -check bounds -traceback -g $(DEFINES) -DFDMDEBUG -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL
else
  FFLAGS = -fpp -I../patch/lagRamses -O3 $(DEFINES) -DLONGINT -DQUADHILBERT -DOUTPUT_PARTICLE_POTENTIAL
endif
ifdef FDMTRACE
  FFLAGS += -DFDMDEBUG
endif
ifeq ($(MG_PROFILE),1)
  # Refined-level multigrid phase timers only; retain the production -O3 flags.
  FFLAGS += -DMG_PROFILE
endif
ifdef AMR_INDEX_CHECK
  FFLAGS += -DAMR_INDEX_CHECK
endif
ifdef MORTON128
  FFLAGS += -DMORTON128
endif
ifeq ($(CIC_ATOMIC_SCATTER),1)
  FFLAGS += -DCIC_ATOMIC_SCATTER
endif

#############################################################################
# HDF5 parallel I/O support (build with: make HDF5=1)
#############################################################################
ifdef HDF5
  HDF5_DIR = $(HOME)/local/hdf5
  FFLAGS += -DHDF5 -I$(HDF5_DIR)/include
  LIBS_HDF5 = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz -Wl,-rpath,$(HDF5_DIR)/lib
  HDF5_MODOBJ = ramses_hdf5_io.o
  HDF5_OBJS = backup_hdf5.o restore_hdf5.o
else
  LIBS_HDF5 =
  HDF5_MODOBJ =
  HDF5_OBJS =
endif

#############################################################################
# CUDA GPU hydro support (build with: make USE_CUDA=1)
#############################################################################
ifdef USE_CUDA
  CUDA_DIR = /opt/ohpc/pub/cuda/13.0.2
  NVCC = $(CUDA_DIR)/bin/nvcc
  CUDA_ARCH = -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90
  N_CUDA_STREAMS ?= 16
  NVCCFLAGS = -O3 $(CUDA_ARCH) -DNDIM=$(NDIM) -DNVAR=$(NVAR) -DNVECTOR=$(NVECTOR) -DN_STREAMS=$(N_CUDA_STREAMS)
  FFLAGS += -DHYDRO_CUDA
  LIBS_CUDA = -L$(CUDA_DIR)/lib64 -lcudart -lcufft -lcublas -lcublasLt -lstdc++ -Wl,-rpath,$(CUDA_DIR)/lib64
  CUDA_MODOBJ = cuda_commons.o hydro_cuda_interface.o poisson_cuda_interface.o sink_cuda_interface.o scalar_cuda_interface.o particle_cuda_interface.o snrt_cuda_interface.o
  CUDA_OBJS = cuda_stream_pool.o hydro_cuda_kernels.o poisson_cuda_kernels.o sink_cuda_kernels.o scalar_cuda_kernels.o particle_cuda_kernels.o snrt_cuda_kernels.o
ifdef USE_CUFFTMP
  CUFFTMP_DIR = $(HOME)/miniconda3/lib/python3.13/site-packages/nvidia/cufftmp/cu12
  FFLAGS += -DUSE_CUFFTMP
  NVCCFLAGS += -DUSE_CUFFTMP
  LIBS_CUDA += -L$(CUFFTMP_DIR)/lib -lcufftMp \
               -Wl,-rpath,$(CUFFTMP_DIR)/lib
endif

else
  LIBS_CUDA =
  CUDA_MODOBJ =
  CUDA_OBJS =
endif

#############################################################################
# Eunha cooling module (always linked; runtime selection via cooling_method)
#############################################################################
EUNHA_MODOBJ = eunha_cooling_mod.o

ifeq ($(SNRT),1)
  ifneq ($(USE_CUDA),1)
    $(error SNRT=1 requires USE_CUDA=1)
  endif
  FFLAGS += -DSNRT
  ifeq ($(DUST_LIVE),1)
    FFLAGS += -DDUST_LIVE
  endif
  # The driver is not self-contained: keep the complete Fortran module graph
  # in the production link. CUDA entry points are supplied by
  # snrt_cuda_kernels.o; these interface modules only publish the C ABI.
  SNRT_MODOBJ = mpi_mod.o snrt_spectral_contract.o snrt_agn_source.o \
                snrt_thermochemistry.o snrt_state.o snrt_agn_locator.o \
                snrt_nlte_coupling.o snrt_angular_quadrature.o \
                snrt_amr_topology.o snrt_transport_step.o \
                snrt_rt_transaction.o snrt_dust_ir.o snrt_dust_coupling.o \
                snrt_dust_transaction.o snrt_dust_receiver.o snrt_dust_contract.o \
                snrt_cuda_ledger_interface.o \
                snrt_cuda_sparse_transport_interface.o \
                snrt_cuda_limited_rt_step_interface.o \
                snrt_cuda_multigroup_interface.o snrt_ramses_driver.o
  ifeq ($(SNRT_LEDGER_DIAGNOSTIC),1)
    FFLAGS += -DSNRT_LEDGER_DIAGNOSTIC
    SNRT_MODOBJ += snrt_agn_ledger.o
  endif
else
  SNRT_MODOBJ =
endif
ifeq ($(DUST_LIVE),1)
  ifneq ($(SNRT),1)
    $(error DUST_LIVE=1 requires SNRT=1)
  endif
endif

#############################################################################
# FFTW3 CPU direct Poisson solver (project default; disable with USE_FFTW=0)
#############################################################################
ifeq ($(USE_FFTW),1)
  FFTW3_DIR = $(HOME)/local
  FFLAGS += -DUSE_FFTW -I$(FFTW3_DIR)/include
  LIBS_FFTW = -L$(FFTW3_DIR)/lib -lfftw3_mpi -lfftw3_omp -lfftw3 -lm \
              -Wl,-rpath,$(FFTW3_DIR)/lib
else
  LIBS_FFTW =
endif

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

LIBGRACKLE = -lgrackle -lhdf5 -lz
LIBS = $(LIBMPI) $(LIBS_HDF5) $(LIBS_CUDA) $(LIBS_FFTW)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../patch/cuda:../patch/oct_tree:../patch/cuRamses:../$(SOLVER):../aton:../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ	= amr_parameters.jaehyun.o amr_commons.kjhan.o timing.o radiation_commons.o amr_index.o random.o pm_parameters.o pm_commons.o poisson_parameters.o \
	 poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o mpi_large.o ksection.o dynamic_exchange.o sparse_mat.o \
	 morton_keys.o morton_hash.o neutrino_commons.o scalar_de_commons.o dark_energy_commons.o power_spectrum.o \
	 clfind_commons.o gadgetreadfile.o write_makefile.o write_gitinfo.o \
	 $(EUNHA_MODOBJ) $(HDF5_MODOBJ) $(CUDA_MODOBJ) $(SNRT_MODOBJ) snrt_agn_efficiency.o agn_feedback_deposition.o dark_cooling_mod.o fdm_commons.o pbh_commons.o \
	 $(STELLAR_ENRICHMENT_MODOBJ)

# SIDM object (depends on amr_commons, pm_commons, random)
SIDM_OBJ = sidm_scatter.o sidm_angular.o sidm_estimator.o
# Atomic Dark Matter objects
ADM_OBJ = dark_cooling_fine.o adm_hpm_mod.o adm_hpm_fine.o
# Evaporating primordial black hole objects
PBH_OBJ = pbh_evap_fine.o
# Fuzzy Dark Matter objects
FDM_OBJ = fdm_step.o fdm_hjm.o output_fdm.o light_cone.fdm.o

AMROBJ = read_params.jaehyun.o init_amr.o init_time.kjhan.o init_refine.o adaptive_loop.jaehyun.o amr_step.jaehyun.o update_time.o \
         output_amr.kjhan.o flag_utils.kjhan.o physical_boundaries.kjhan.o virtual_boundaries.kjhan.o refine_utils.o nbors_utils.kjhan.o \
         hilbert.o load_balance.kjhan.o title.o sort.o cooling_fine.kjhan.o sgs_fine.o diag_eint.o units.o light_cone.part.o light_cone.hydro2.o light_cone.sink2.o movie_mod.yonghwi_org.o kjhan.o\
		 output_sphere_hydro.o output_sphere_part.o output_sphere_sink.o \
         morton_init.o $(HDF5_OBJS)
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.kjhan.o synchro_fine.kjhan.o move_fine.kjhan.o newdt_fine.kjhan.o particle_tree.kjhan.o \
        add_list.o remove_list.o star_formation.kjhan.o sink_particle.kjhan.o stellar_ramses_runtime.o feedback.kjhan3.o \
		clump_finder.o clump_merger.o \
        flag_formation_sites.o init_sink.o output_sink.o
# Poisson solver objects
#POISSONOBJ = init_poisson.o phi_fine_cg.kisti.o interpol_phi.kjhan.o force_fine.kjhan.o multigrid_coarse.kjhan.o multigrid_fine_commons.o \
#             multigrid_fine_fine.kisti.o multigrid_fine_coarse.kisti.o gravana.o boundary_potential.o rho_ana.kjhan.o output_poisson.o
POISSONOBJ = init_poisson.o phi_fine_cg.kjhan.o interpol_phi.kjhan.o force_fine.kjhan.o multigrid_coarse.kjhan.o multigrid_fine_commons.o \
             multigrid_fine_fine.kjhan.o multigrid_fine_coarse.kjhan.o gravana.o boundary_potential.kjhan.o rho_ana.kjhan.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.kjhan.o godunov_fine.kjhan.o \
           uplmde.kjhan.o umuscl.kjhan.o interpol_hydro.kjhan.o godunov_utils.kjhan.o condinit.o hydro_flag.kjhan.o hydro_boundary.o \
           boundana.o read_hydro_params.o synchro_hydro_fine.kjhan.o
# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) $(SIDM_OBJ) $(ADM_OBJ) $(PBH_OBJ) $(FDM_OBJ)
# ATON objects
ATON_MODOBJ = rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o rad_backup.o ../aton/atonlib/libaton.a

# Phase 0 Fortran-native stellar enrichment modules.
PHASE0_STELLAR_ENRICHMENT ?= 1
ifeq ($(PHASE0_STELLAR_ENRICHMENT),1)
  FFLAGS += -DPHASE0_STELLAR_ENRICHMENT
endif
STELLAR_YIELD_TABLE ?=
STELLAR_EMBEDDED_SOURCE = $(PATCH)/stellar_yield_embedded_data.f90
ifdef STELLAR_EMBEDDED_YIELDS
  FFLAGS += -DSTELLAR_EMBEDDED_YIELDS
  STELLAR_EMBEDDED_MODOBJ = stellar_yield_embedded_data.o
else
  STELLAR_EMBEDDED_MODOBJ =
endif
# The production source-order sidecar is load-bearing.  The RAMSES runtime
# module is listed last in the sidecar but remains in PMOBJ below so the
# integration call site and its module dependencies stay in the production
# link order.
include $(PATCH)/stellar_enrichment_sources.mk
PHASE0_STELLAR_ENRICHMENT_OBJECTS := $(notdir $(PHASE0_STELLAR_ENRICHMENT_SOURCES:.f90=.o))
STELLAR_ENRICHMENT_MODOBJ = $(filter-out stellar_ramses_runtime.o,$(PHASE0_STELLAR_ENRICHMENT_OBJECTS))
ifdef STELLAR_EMBEDDED_YIELDS
  STELLAR_ENRICHMENT_MODOBJ += $(STELLAR_EMBEDDED_MODOBJ)
endif
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) $(CUDA_OBJS) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) $(CUDA_OBJS) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm -f $(WRITE_MAKEFILE_SOURCE)
#	rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	rm -f $(WRITE_MAKEFILE_SOURCE)
#	rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
 -DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@
WRITE_MAKEFILE_SOURCE = write_makefile.generated.f90
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(firstword $(MAKEFILE_LIST)) $(WRITE_MAKEFILE_SOURCE)
	$(F90) $(FFLAGS) -c $(WRITE_MAKEFILE_SOURCE) -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $@
# ifx ICE workaround: phi_fine_cg.kisti crashes with -qopenmp
phi_fine_cg.kisti.o: phi_fine_cg.kisti.f90
	mpiifx $(FFLAGS) -c $^ -o $@
ifdef USE_CUDA
%.o:%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
hydro_cuda_kernels.o: cuda_stream_pool.o
poisson_cuda_kernels.o: cuda_stream_pool.o ../patch/cuRamses/amr_cuda_index.cuh
sink_cuda_kernels.o: cuda_stream_pool.o
scalar_cuda_kernels.o: cuda_stream_pool.o ../patch/cuRamses/amr_cuda_index.cuh
particle_cuda_kernels.o: cuda_stream_pool.o ../patch/cuRamses/amr_cuda_index.cuh
hydro_cuda_interface.o: cuda_commons.o
poisson_cuda_interface.o: cuda_commons.o
sink_cuda_interface.o: cuda_commons.o
scalar_cuda_interface.o: cuda_commons.o
particle_cuda_interface.o: cuda_commons.o
cuda_commons.o: amr_commons.kjhan.o
rho_fine.kjhan.o move_fine.kjhan.o synchro_fine.kjhan.o: particle_cuda_interface.o
force_fine.kjhan.o: scalar_cuda_interface.o
multigrid_fine_commons.o: poisson_cuda_interface.o
endif
amr_index.o: amr_commons.kjhan.o
# SNRT module graph.  Explicit prerequisites are required because the
# implicit %.o rule otherwise tries to compile a module before its .mod file
# exists when a clean production build is requested.
snrt_spectral_contract.o: amr_parameters.jaehyun.o
snrt_state.o: amr_commons.kjhan.o snrt_spectral_contract.o snrt_thermochemistry.o
snrt_agn_locator.o: amr_commons.kjhan.o
snrt_agn_source.o: amr_parameters.jaehyun.o
snrt_nlte_coupling.o: snrt_agn_source.o
snrt_thermochemistry.o: amr_parameters.jaehyun.o snrt_agn_source.o
snrt_angular_quadrature.o: snrt_state.o
snrt_amr_topology.o: amr_commons.kjhan.o mpi_mod.o snrt_state.o
snrt_transport_step.o: snrt_state.o snrt_angular_quadrature.o \
                        snrt_amr_topology.o \
                        mpi_mod.o \
                        snrt_cuda_sparse_transport_interface.o \
                        snrt_cuda_limited_rt_step_interface.o \
                        snrt_cuda_multigroup_interface.o
snrt_rt_transaction.o: amr_parameters.jaehyun.o mpi_mod.o
snrt_dust_ir.o: ../patch/lagRamses/snrt_dust_ir.f90
snrt_dust_coupling.o: snrt_thermochemistry.o ../patch/lagRamses/snrt_dust_coupling.f90
snrt_dust_transaction.o: amr_parameters.jaehyun.o
snrt_dust_receiver.o: ../patch/lagRamses/snrt_dust_receiver.f90
snrt_dust_contract.o: ../patch/lagRamses/snrt_dust_contract.f90
snrt_agn_efficiency.o: amr_parameters.jaehyun.o
agn_feedback_deposition.o: amr_parameters.jaehyun.o
init_sink.o: pm_commons.o
sink_particle.kjhan.o: snrt_agn_efficiency.o agn_feedback_deposition.o amr_parameters.jaehyun.o pm_commons.o
ifeq ($(PHASE0_STELLAR_ENRICHMENT),1)
sink_particle.kjhan.o: stellar_enrichment_config.o
endif
read_params.jaehyun.o: snrt_agn_efficiency.o
read_params.jaehyun.o: snrt_thermochemistry.o snrt_dust_contract.o
ifeq ($(SNRT),1)
read_params.jaehyun.o: snrt_spectral_contract.o
endif
snrt_ramses_driver.o: amr_commons.kjhan.o hydro_parameters.o hydro_commons.o \
                       pm_commons.o snrt_agn_efficiency.o \
                       snrt_state.o snrt_amr_topology.o snrt_transport_step.o \
                       snrt_rt_transaction.o \
                       snrt_dust_transaction.o snrt_dust_receiver.o snrt_dust_contract.o \
                       snrt_angular_quadrature.o snrt_agn_locator.o \
                       snrt_agn_source.o snrt_nlte_coupling.o snrt_thermochemistry.o \
                       snrt_spectral_contract.o \
                       snrt_cuda_interface.o snrt_cuda_ledger_interface.o \
                       snrt_cuda_sparse_transport_interface.o \
                       snrt_cuda_limited_rt_step_interface.o \
                       snrt_cuda_multigroup_interface.o
snrt_agn_ledger.o: amr_commons.kjhan.o pm_commons.o
# [RESIZABLE] grow_grid_capacity preserves every capacity bundle, including
# the optional radiation fields; build those modules before pm_commons.mod.
pm_commons.o: amr_commons.kjhan.o pm_parameters.o hydro_commons.o poisson_commons.o morton_hash.o radiation_commons.o
radiation_commons.o: amr_commons.kjhan.o timing.o
morton_hash.o: morton_keys.o amr_index.o
hydro_commons.o: hydro_parameters.o
poisson_commons.o: poisson_parameters.o
stellar_enrichment_contract.o: stellar_enrichment_config.o
stellar_snia_dtd.o: stellar_enrichment_config.o
stellar_snia_population_contract.o: stellar_enrichment_config.o stellar_snia_dtd.o
stellar_snia_physical_contract.o: stellar_enrichment_config.o
stellar_snia_event_ledger.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_snia_physical_contract.o
stellar_native_units.o: stellar_enrichment_config.o
stellar_snia_cell_deposition.o: stellar_enrichment_config.o stellar_snia_physical_contract.o stellar_native_units.o
stellar_progress_contract.o: stellar_enrichment_config.o
stellar_snia_runtime_accounting.o: stellar_enrichment_config.o
stellar_yield_tables.o: stellar_enrichment_config.o
stellar_yield_backend.o: stellar_yield_tables.o $(STELLAR_EMBEDDED_MODOBJ)
stellar_yield_embedded_data.o: $(STELLAR_EMBEDDED_SOURCE) stellar_yield_tables.o stellar_enrichment_config.o
stellar_yield_interpolation.o: stellar_enrichment_config.o stellar_yield_tables.o
stellar_yield_provider.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_yield_tables.o stellar_yield_interpolation.o
stellar_ssp_sources.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_yield_tables.o stellar_yield_provider.o
stellar_source_increment.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_yield_tables.o stellar_ssp_sources.o
stellar_population_ledger.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_ssp_sources.o stellar_snia_physical_contract.o
stellar_enrichment_driver.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_yield_tables.o stellar_ssp_sources.o stellar_source_increment.o stellar_population_ledger.o
stellar_cell_deposition.o: stellar_enrichment_config.o stellar_enrichment_contract.o
stellar_ramses_field_map.o: stellar_enrichment_config.o
stellar_ramses_bridge.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_snia_physical_contract.o stellar_snia_cell_deposition.o stellar_ramses_field_map.o
stellar_ramses_mapped_bridge.o: stellar_enrichment_config.o stellar_enrichment_contract.o stellar_ramses_field_map.o
stellar_ramses_runtime.o: amr_commons.kjhan.o pm_commons.o hydro_commons.o stellar_enrichment_config.o stellar_enrichment_contract.o stellar_native_units.o stellar_progress_contract.o stellar_snia_runtime_accounting.o stellar_yield_tables.o stellar_yield_backend.o stellar_yield_audit.o stellar_enrichment_driver.o stellar_population_ledger.o stellar_snia_population_contract.o stellar_snia_physical_contract.o stellar_snia_cell_deposition.o stellar_ramses_field_map.o stellar_ramses_bridge.o
feedback.kjhan3.o: stellar_ramses_runtime.o
read_hydro_params.o: stellar_enrichment_config.o
output_amr.kjhan.o: stellar_enrichment_config.o stellar_ramses_runtime.o
ramses.o: stellar_enrichment_config.o stellar_ramses_runtime.o
amr_step.jaehyun.o: stellar_enrichment_config.o
stellar_yield_audit.o: stellar_enrichment_config.o stellar_yield_tables.o

$(STELLAR_EMBEDDED_SOURCE):
	@test -n "$(STELLAR_YIELD_TABLE)" || (echo "Set STELLAR_YIELD_TABLE to generate the embedded yield module"; exit 2)
	python3 $(PATCH)/generate_stellar_yield_module.py --input "$(STELLAR_YIELD_TABLE)" --output "$@"
power_spectrum.o: amr_parameters.jaehyun.o amr_commons.kjhan.o poisson_commons.o
scalar_de_commons.o: amr_parameters.jaehyun.o amr_commons.kjhan.o
dark_energy_commons.o: scalar_de_commons.o
init_time.o: scalar_de_commons.o
multigrid_fine_commons.o: scalar_de_commons.o dark_energy_commons.o
synchro_fine.kjhan.o: scalar_de_commons.o
output_fdm.o: amr_commons.kjhan.o poisson_commons.o ksection.o morton_keys.o morton_hash.o fdm_commons.o
output_amr.kjhan.o: amr_parameters.jaehyun.o amr_commons.kjhan.o
dynamic_exchange.o: amr_parameters.jaehyun.o amr_commons.kjhan.o mpi_large.o
init_part.o: dynamic_exchange.o
init_hydro.o: dynamic_exchange.o
virtual_boundaries.kjhan.o: dynamic_exchange.o
dynamic_exchange_mpi_test: amr_parameters.jaehyun.o amr_commons.kjhan.o pm_parameters.o pm_commons.o bisection.o mpi_large.o ksection.o dynamic_exchange.o ../tests/dynamic_exchange_mpi.f90
	$(F90) $(FFLAGS) amr_parameters.jaehyun.o amr_commons.kjhan.o pm_parameters.o pm_commons.o bisection.o mpi_large.o ksection.o dynamic_exchange.o \
		../tests/dynamic_exchange_mpi.f90 -o $@
light_cone.fdm.o: amr_commons.kjhan.o pm_commons.o poisson_commons.o
eunha_cooling_mod.o: amr_parameters.jaehyun.o amr_commons.kjhan.o
dark_cooling_mod.o: amr_parameters.jaehyun.o
dark_cooling_fine.o: amr_commons.kjhan.o pm_commons.o dark_cooling_mod.o
FORCE:

.PHONY: snrt_bundle_gate
snrt_bundle_gate:
	../simulation/snrt/tests/run_snrt_bundle_gate.sh

#############################################################################
clean :
	rm -f *.o *.$(MOD)
#############################################################################
