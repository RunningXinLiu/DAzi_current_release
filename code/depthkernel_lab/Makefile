FC = gfortran
FFLAGS = -O2 -fopenmp -ffixed-line-length-none -ffree-line-length-none -ffloat-store -fbounds-check -m64

SRC_DIR = ../src_inv_iso_joint
PROD_SRC_DIR = ../src_inv_iso_joint
PYTHON_DISBA = /opt/miniconda3/envs/dispa/bin/python
NUMBA_CACHE_DIR = /private/tmp/numba_cache

all: depthkernel_baseline depthkernel_prod1d_driver depthkernel_ti_worker surfdisp96_forward_batch

depthkernel_baseline: depthkernel_baseline.o surfdisp96.o
	$(FC) -fopenmp $^ -o $@

depthkernel_baseline.o: depthkernel_baseline.f90
	$(FC) $(FFLAGS) -c $< -o $@

surfdisp96.o: $(SRC_DIR)/surfdisp96.f
	$(FC) $(FFLAGS) -c $< -o $@

depthkernel_prod1d_driver: depthkernel_prod1d_driver.o DepthKernel1D.o prod_surfdisp96.o
	$(FC) -fopenmp $^ -o $@

depthkernel_ti_worker: depthkernel_ti_worker.o prod_surfdisp96.o prod_tregn96.o
	$(FC) $^ -o $@

surfdisp96_forward_batch: surfdisp96_forward_batch.o prod_surfdisp96.o
	$(FC) -fopenmp $^ -o $@

depthkernel_prod1d_driver.o: depthkernel_prod1d_driver.f90
	$(FC) $(FFLAGS) -c $< -o $@

depthkernel_ti_worker.o: depthkernel_ti_worker.f90
	$(FC) $(FFLAGS) -c $< -o $@

surfdisp96_forward_batch.o: surfdisp96_forward_batch.f90
	$(FC) $(FFLAGS) -c $< -o $@

DepthKernel1D.o: $(PROD_SRC_DIR)/DepthKernel1D.f90
	$(FC) $(FFLAGS) -c $< -o $@

prod_surfdisp96.o: $(PROD_SRC_DIR)/surfdisp96.f
	$(FC) $(FFLAGS) -c $< -o $@

prod_tregn96.o: $(PROD_SRC_DIR)/tregn96.f
	$(FC) $(FFLAGS) -c $< -o $@

run: depthkernel_baseline
	./depthkernel_baseline sample_column.in baseline_kernel.out

prod1d: depthkernel_prod1d_driver
	./depthkernel_prod1d_driver sample_column.in prod1d_kernel.out

prod-compare: run prod1d
	grep -v '^# depthkernel_' baseline_kernel.out | grep -v '^# cpu_seconds' > baseline_kernel.norm
	grep -v '^# depthkernel_' prod1d_kernel.out | grep -v '^# cpu_seconds' > prod1d_kernel.norm
	diff -u baseline_kernel.norm prod1d_kernel.norm > prod1d_vs_baseline.diff || test $$? -eq 1
	@if [ -s prod1d_vs_baseline.diff ]; then \
		echo "Production 1D kernel differs from baseline; see prod1d_vs_baseline.diff"; \
		exit 1; \
	else \
		echo "Production 1D kernel matches baseline values."; \
	fi

prod-real-compare: depthkernel_baseline depthkernel_prod1d_driver
	$(PYTHON_DISBA) compare_prod1d_real_columns.py

prod-real-modref-compare: depthkernel_baseline depthkernel_prod1d_driver
	$(PYTHON_DISBA) compare_prod1d_real_columns.py --mod /Users/liuxin/Desktop/DSurf_test/DAzi_large_data_260128/example/anning_real_SRL_ani_cell30_15_nreal050_iaratio0p5/MOD_Ref --outdir real_column_prod1d_regression_modref --dln-values 0.02,0.01,0.005 --minthk-values 1,2,4

compare: run
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) compare_disba.py --input sample_column.in --fortran baseline_kernel.out --compare-out disba_compare.out --layer-out disba_layer_kernels.out

batch-compare: run
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) batch_disba_fd.py --input sample_column.in --fortran baseline_kernel.out --compare-out batch_disba_compare.out --workers 4

pdsurf-numba-compare: run
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) batch_pdsurf_numba_fd.py --input sample_column.in --fortran baseline_kernel.out --compare-out pdsurf_numba_compare.out

pdsurf-tile-benchmark:
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) benchmark_pdsurf_tile_batch.py

sweep: depthkernel_baseline
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) sweep_kernel_params.py

sweep-real: depthkernel_baseline
	NUMBA_CACHE_DIR=$(NUMBA_CACHE_DIR) $(PYTHON_DISBA) sweep_real_columns.py

summarize:
	$(PYTHON_DISBA) summarize_sweep.py --summary real_column_sweep_modref/real_columns_sweep_summary.csv --outdir real_column_sweep_modref/summary_plots

clean:
	rm -f *.o *.mod depthkernel_baseline depthkernel_prod1d_driver depthkernel_ti_worker surfdisp96_forward_batch baseline_kernel.out prod1d_kernel.out baseline_kernel.norm prod1d_kernel.norm prod1d_vs_baseline.diff disba_compare.out disba_layer_kernels.out sweep_summary.csv
