# Phase 0 Fortran source order.
#
# Include this list from the lagRamses build configuration.  The order is the
# module dependency order and keeps all new stellar-enrichment source in this
# patch directory.

PHASE0_STELLAR_ENRICHMENT_SOURCES = \
  patch/lagRamses/stellar_enrichment_config.f90 \
  patch/lagRamses/stellar_enrichment_contract.f90 \
  patch/lagRamses/stellar_yield_tables.f90 \
  patch/lagRamses/stellar_yield_interpolation.f90 \
  patch/lagRamses/stellar_yield_provider.f90 \
  patch/lagRamses/stellar_ssp_sources.f90 \
  patch/lagRamses/stellar_source_increment.f90 \
  patch/lagRamses/stellar_enrichment_driver.f90 \
  patch/lagRamses/stellar_cell_deposition.f90 \
  patch/lagRamses/stellar_ramses_field_map.f90 \
  patch/lagRamses/stellar_ramses_bridge.f90 \
  patch/lagRamses/stellar_ramses_mapped_bridge.f90 \
  patch/lagRamses/stellar_ramses_runtime.f90 \
  patch/lagRamses/stellar_yield_audit.f90
