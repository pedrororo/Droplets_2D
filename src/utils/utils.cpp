#include "utils.h"

void usage (const char pname[])
{
	printf("%s\n",PRINT_LINE);
	printf("Usage:> %s <configuration_file>\n",pname);
	printf("%s\n",PRINT_LINE);
	printf("\t<configuration_file> = Configuration file with the parameter values from the network\n");
	printf("%s\n",PRINT_LINE);
	printf("Example:\n");
	printf("\t%s example_configs/create_spiral.ini\n",pname);
	printf("%s\n",PRINT_LINE);
}
