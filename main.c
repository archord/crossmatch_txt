/*
######################################################################
## Time-stamp: 
## Filename:      $Name:  $       
## Version:       $Revision: 1.2 $ 
## Author:        Yang Xu <yxuctgu@gmail.com>
## Purpose:       cross match of astrometry.
## CVSAuthor:     $Author: cyxu $ 
## Note:          
#-                
## $Id: main.cpp,v 1.2 2012/04/16 08:11:15 cyxu Exp $
#======================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "function.h"

extern int showResult;		//0输出所有结果，1输出匹配的结果，2输出不匹配的结果
extern int printResult;		//0将匹配结果不输出到终端，1输出到终端
extern int showProcessInfo;	//0不输出处理过程信息，1输出
extern int fitsHDU;         //fitsHDU=3: read fits file from the 3rd hdu
extern int useCross;        //use cross method
extern double minZoneLength;//the min length of the zone's side
extern double searchRadius;
extern double areaBox;      //判断两颗星是一颗星的最大误差，默认为20角秒

extern int dbConfigInCommandLine;
extern char *configFile;
extern char *host;
extern char *port;
extern char *dbname;
extern char *user;
extern char *password;

extern char *catfile_table;
extern char *match_table;
extern char *ot_table;
extern char *ot_flux_table;

extern int areaWidth;
extern int areaHeight;
extern float planeErrorRedius;

extern int fluxRatioSDTimes; //factor of flux filter

void showHelp(){
	printf("usage: crossmatch refFile sampleFile outFile errorBox\n");
}


int main(int argc, char** argv){

	if(argc != 5){
		showHelp();
		return 0;
	}
    
    areaBox = atof(argv[4]);

    mainSphere(argv[1],argv[2],argv[3]);
    
	return 0;
}
