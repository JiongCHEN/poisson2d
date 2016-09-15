EXE=../build/examples/main
XDIM=16
YDIM=16
OUTDIR=../result/grid-$(XDIM)-$(YDIM)

$(OUTDIR)/grid.dat: $(EXE)
	mkdir -p $(OUTDIR)
	$(EXE) $(OUTDIR)/grid.dat
