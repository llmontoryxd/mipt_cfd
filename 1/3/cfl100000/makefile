CFLAGS = -fpconstant -openmp -O3 -r8 -fpp -D_metis #  -DMpiMultiBlock  
COMPILER = mpiifort

objects = core1 core2 flowcore flowio m_perfect flowsolver driver
                 
result: $(objects)
	$(COMPILER) $(CFLAGS)  -o  flowmod2017mpi *.o  -L./ -lmetis -lm 

core1:  NesvetayCore/bc_types.f90 NesvetayCore/common_proc.f90  NesvetayCore/basic_classes.f90 NesvetayCore/reconst_lsq.f90 NesvetayCore/reconst_tvd1d.f90
	$(COMPILER) $(CFLAGS)  -c  $^

core2: NesvetayCore/volume_mesh_class.f90 NesvetayCore/mesh_reading.f90  NesvetayCore/mpi_ghosts.f90  NesvetayCore/mesher_io.f90 NesvetayCore/nesvetay_drv.f90
	$(COMPILER) $(CFLAGS)  -c  $^

flowcore:  FlowCore/flow_class.f90 FlowCore/flow_common.f90  FlowCore/rs_common.f90  FlowCore/numerics_common.f90 FlowCore/bc_common.f90 FlowCore/restart.f90
	$(COMPILER) $(CFLAGS)  -c  $^

flowio:   FlowIO/revolve_io.f90 FlowIO/line_io.f90 FlowIO/surface_io.f90  FlowIO/volume_io.f90  FlowIO/vtk_io.f90 FlowIO/flow_io.f90
	$(COMPILER) $(CFLAGS)  -c  $^

m_perfect:  Models/perfect/model_perfect.f90 Models/perfect/rs_perfect.f90 Models/perfect/numerics_perfect.f90 Models/perfect/bc_perfect.f90 Models/perfect/init_perfect.f90 
	$(COMPILER) $(CFLAGS)  -c  $^

flowsolver: FlowSolver/partitioning.f90  FlowSolver/prim_reconst.f90 FlowSolver/scheme_rhs.f90 FlowSolver/initialize.f90 FlowSolver/lu_sgs.f90 FlowSolver/lu_sgs_multi.f90
	$(COMPILER) $(CFLAGS)  -c  $^


driver:  Driver/steady_solver.f90   Driver/driver.f90 
	$(COMPILER) $(CFLAGS)  -c  $^

clean:
	rm *.o
	rm *.mod
