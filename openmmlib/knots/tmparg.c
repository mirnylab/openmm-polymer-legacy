
static char *config_file = NULL;
static char *input_file  = NULL;

#define USAGE "-f <filein> -c <configfile>"
#define ARGS "f:c:"

int parse_args(int argc, char **argv){
  int c,errflg = 0;

  optarg = NULL;
  optind = 0;

  while(!errflg && (c=getopt(argc,argv,ARGS))!=-1){
    switch(c){
      case 'c':
	config_file = optarg;
	break;
      case 'f':
	input_file  = optarg;
	break;
    }
  }

  if(optind < 1 || argc < 2)
    errflg=1;

  if(errflg){
    fprintf(stderr,"%s\nUsage: %s %s\n", DISCLAIMER, argv[0], USAGE);
    exit(0);
  }

  filein = argv[optind];

  return 0;
}




