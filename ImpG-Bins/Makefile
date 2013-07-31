CC = g++
COPTS = -O3 -lm -Wall

BINS = ImpG-Summary-GenBeta ImpG-Summary \
	GenLD ImpG-SummaryLD-GenBeta ImpG-SummaryLD #\
#	ImpG-Raw-GenBeta

all: $(BINS)

# ImpG-Raw-GenBeta
#ImpG-Raw-GenBeta: impg.o impg_raw_gen_beta_main.o util.o linsubs.o
#	$(CC) impg.o impg_raw_gen_beta_main.o util.o linsubs.o \
#		-o ImpG-Raw-GenBeta $(COPTS)
#impg_raw_gen_beta_main.o: impg_raw_gen_beta_main.cpp
#	$(CC) -c impg_raw_gen_beta_main.cpp $(COPTS)

#GenLD
GenLD: genld_main.o genld.o util.o
	$(CC) genld_main.o genld.o util.o -o GenLD $(COPTS)
genld_main.o: genld_main.cpp
	$(CC) -c genld_main.cpp $(COPTS)
genld.o: genld.cpp genld.h
	$(CC) -c genld.cpp $(COPTS)

# ImpG-SummaryLD
ImpG-SummaryLD: impg.o impz_ld_main.o util.o linsubs.o
	$(CC) impg.o impz_ld_main.o util.o linsubs.o \
		-o ImpG-SummaryLD $(COPTS)
impz_ld_main.o: impz_ld_main.cpp
	$(CC) -c impz_ld_main.cpp $(COPTS)
	
# ImpG-Summary
ImpG-Summary: impg.o impz_main.o util.o linsubs.o
	$(CC) impg.o impz_main.o util.o linsubs.o \
		-o ImpG-Summary $(COPTS)
impg_summary_main.o: impz_main.cpp
	$(CC) -c impz_main.cpp $(COPTS)

# ImpG-SummaryLD-GenBeta
ImpG-SummaryLD-GenBeta: impg.o zgenbt_ld_main.o util.o linsubs.o
	$(CC) impg.o zgenbt_ld_main.o util.o linsubs.o \
		-o ImpG-SummaryLD-GenBeta $(COPTS)
zgenbt_ld_main.o: zgenbt_ld_main.cpp
	$(CC) -c zgenbt_ld_main.cpp $(COPTS)
	
# ImpG-Summary-GenBeta
ImpG-Summary-GenBeta: impg.o zgenbt_main.o util.o linsubs.o
	$(CC) impg.o zgenbt_main.o util.o linsubs.o \
		-o ImpG-Summary-GenBeta $(COPTS)
zgenbt_main.o: zgenbt_main.cpp
	$(CC) -c zgenbt_main.cpp $(COPTS)
	
# ImpG-Summary commons
impg.o: impg.cpp impg.h
	$(CC) -c impg.cpp $(COPTS)

# all common
util.o: util.cpp util.h
	$(CC) -c util.cpp $(COPTS)
linsubs.o: linsubs.cpp linsubs.h
	$(CC) -c linsubs.cpp $(COPTS)

clean:
	rm -f *.o $(BINS)
