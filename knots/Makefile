CKNOT=CODE_KNOT_ANALYSE

SCRIPTS=knotsubmit.pl getreduced.pl getcif.pl analyze_structure_kgb.pl KGBCIF.pm 
DWNSCRIPT=knotsubmit.pl

BINARIES=$(CKNOT)/Reduce_knot20

CGIDIR=/var/www/cgi-bin/knots
HTMLDIR=/var/www/html/knots

CFLAGS=-m32 -Ofast -m32  

HFILES=html/index.html html/3_1.jpg html/4_1.jpg html/5_1.jpg html/5_2.jpg  html/question.png html/knot_full_small.jpg html/knot_small.jpg

all: htmlstuff cknot


web: cknot htmlstuff
	/bin/cp -f $(SCRIPTS) $(BINARIES) $(CGIDIR)
	/bin/cp -f $(HFILES) $(HTMLDIR)
	cd $(CGIDIR) && ln -f $(DWNSCRIPT) knotsout.zip

cknot:
	cd $(CKNOT) && make

htmlstuff:
	cd html && make
