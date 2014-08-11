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
## $Id: function.cpp,v 1.2 2012/04/16 08:11:15 cyxu Exp $
#======================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "function.h"

#define MaxStringLength 255

extern struct WorldCoor *GetFITSWCS(char *filename,
        char *header,
        int verbose,
        double *cra,
        double *cdec,
        double *dra,
        double *ddec,
        double *secpix,
        int *wp,
        int *hp,
        int *sysout,
        double *eqout);

float getAngleFromGreatCircle(double dec, double errorRadius);

/*
 *算法说明：
 *1，按照天文坐标，赤经（0~360），赤纬（-90~+90），按每一度为单位将整个球面坐标分成360*180区域，每个区域的大小为赤经1度，赤纬1度。
 *2，建立一个大小为360*180的AREANODE结构体数组（后面称该数组为树干），按照赤经为行，赤纬为列的行列顺序进行存储，即180行，360列，
 **  每个AREANODE结构体中有一个指针，该指针指向该结构体所对应区域中的所有点所组成的链表（后面称该数组为树枝）。
 *3，从源文件中读取数据点，每个数据点对应一个SAMPLE结构体，将该数据点添加到对应的区域链表当中（树干中的对应树枝），
 **  区域链表（树枝）按平衡四叉树进行排列。
 **  例如：点（95.25988，61.11838），对应树干（共360*180个区域）中的第360*61+95个区域（赤经95，赤纬61），将该点添加到该区域的子链表
 **  （树枝）中，该子链表（树枝）按平衡四叉树进行排列。
 *4，从目标文件中读出需要匹配的数据点，先找到该点在树干中对应的区域，再在该区域的树枝中找到对应匹配的点。
 *5，将匹配成功的点对添加到结果链表中去。
 *6，将结果链表输出到文件。
 */

int showResult = 0; //0输出所有结果，1输出匹配的结果，2输出不匹配的结果
int printResult = 0; //0将匹配结果不输出到终端，1输出到终端
int showProcessInfo = 0; //0不输出处理过程信息，1输出
int fitsHDU = 3; //fitsHDU=3: read fits file from the 3rd hdu
int useCross = 0;
double areaBox = ERROR_GREAT_CIRCLE; //判断两颗星是一颗星的最大误差，默认为20角秒
double minZoneLength = ERROR_GREAT_CIRCLE; //3 times of areaBox
double searchRadius = ERROR_GREAT_CIRCLE; //search radius, great circle

/*database configure information*/
int dbConfigInCommandLine = 0; //use command line config databaseinfo or use file
char *configFile = NULL;
char *host = NULL;
char *port = NULL;
char *dbname = NULL;
char *user = NULL;
char *password = NULL;
char *options = NULL;
char *tty = NULL;

/*table information*/
char *catfile_table = NULL;
char *match_table = NULL;
char *ot_table = NULL;
char *ot_flux_table = NULL;

double airmass = 0.0;
double jd = 0.0;

//for sphere coordiante's partition
long raMini;
long decMini;
long raMaxi;
long decMaxi;

float raMinf;
float decMinf;
float raMaxf;
float decMaxf;

int absDecMin; //in north or south, the max is different,
int absDecMax;

int decNode; //number of subarea in dec
int raNode; //number of subarea in ra
int zoneLength; //arc second, length of subarea's width or height
double factor; //=3600/zoneLength

float *raRadiusIndex = NULL;

//for plane coordiante's partition
int areaWidth = 0;
int areaHeight = 0;
float zoneInterval = 0;
int planeZoneX = 0;
int planeZoneY = 0;
float planeErrorRedius = 0.0;

int colNum = 0;

//to filter the flus
double standardDeviation = 0.0;
int fluxRatioSDTimes = 0; //abs(ratioRatio - flusRatioAverage) > flusRatioSD * standardDeviation;
double fluxRatioAverage = 0.0;
double fluxRatioMedian = 0.0;



void printAreaInfo() {
    printf("min ra : %f\n", raMinf);
    printf("min dec: %f\n", decMinf);
    printf("max ra : %f\n", raMaxf);
    printf("max dec: %f\n", decMaxf);
    printf("zone subarea length: %d arc second\n", zoneLength);
    printf("subarea columns(ra): %d\n", raNode);
    printf("subarea rows(dec)  : %d\n", decNode);
}

//init search index radius for the sample's match. 
//though each point's search radius are the same, 
//but the search radius in ra direction is different. 
//we should convert great circle distance to ra angle.
//see function getAngleFromGreatCircle for more detail.

void initRaRadiusIndex() {

    if (decMini > 0 && decMaxi > 0) {
        absDecMin = decMini;
        absDecMax = decMaxi;
    } else if (decMini < 0 && decMaxi < 0) {
        absDecMin = abs(decMaxi);
        absDecMax = abs(decMini);
    } else { // if(decMini<0 && decMaxi>0)
        absDecMin = 0;
        if (abs(raMini) < abs(decMaxi)) {
            absDecMax = abs(decMaxi);
        } else {
            absDecMax = abs(decMini);
        }
    }

    long num = ceil((absDecMax - absDecMin) / searchRadius);
    raRadiusIndex = (float *) malloc(num * sizeof (long));
    long i = 0;
    float tmpDec = 0.0;
    for (i = 0; i < num; i++) {
        tmpDec = absDecMin + i*searchRadius;
        raRadiusIndex[i] = getAngleFromGreatCircle(tmpDec, searchRadius);
    }

    if (showProcessInfo) {
        printf("ra radius index length: %d\n", num);
    }
}

void getAreaBoundary(struct SAMPLE *head) {

    struct SAMPLE *listHead = head;
    struct SAMPLE *tmp = listHead->next; //一般第一个节点为头结点，不存储数据，第二个节点才是第一个数据点
    float raMin;
    float decMin;
    float raMax;
    float decMax;
    if (tmp != NULL) {
        raMin = tmp->alpha;
        raMax = tmp->alpha;
        decMin = tmp->delta;
        decMax = tmp->delta;
        tmp = tmp->next;
    }
    while (tmp) {
        if (tmp->alpha > raMax) {
            raMax = tmp->alpha;
        } else if (tmp->alpha < raMin) {
            raMin = tmp->alpha;
        }
        if (tmp->delta > decMax) {
            decMax = tmp->delta;
        } else if (tmp->delta < decMin) {
            decMin = tmp->delta;
        }
        tmp = tmp->next;
    }

    raMinf = raMin;
    decMinf = decMin;
    raMaxf = raMax;
    decMaxf = decMax;

    decMini = floor(decMin - areaBox);
    decMaxi = ceil(decMax + areaBox);

    float maxDec = fabs(decMaxi) > fabs(decMini) ? fabs(decMaxi) : fabs(decMini);
    float areaBoxForRa = getAngleFromGreatCircle(maxDec, areaBox);
    raMini = floor(raMin - areaBoxForRa);
    raMaxi = ceil(raMax + areaBoxForRa);
}

//get zone's width and hight, width=hight

void getZoneLength() {

    long totalNode = (INDEX_SIZE) / sizeof (struct AREANODE);
    float zoneLengthf = sqrt((decMaxi - decMini + 1)*(raMaxi - raMini + 1)*3600.0 * 3600.0 / totalNode);
    zoneLength = ceil(zoneLengthf);
    if (zoneLength < minZoneLength * 3600)
        zoneLength = ceil(minZoneLength * 3600);
    factor = 3600 / zoneLength;

    decNode = ceil((decMaxi - decMini + 1) * factor);
    raNode = ceil((raMaxi - raMini + 1) * factor);
}

/********************
 *功能：初始化分区信息
 *输入：
 *输出：
 */
struct AREANODE *initAreaNode(struct SAMPLE *point) {

    getAreaBoundary(point);
    getZoneLength();
    initRaRadiusIndex();

    int totalNode = decNode*raNode;
    struct AREANODE *areaTree = (struct AREANODE *) malloc(sizeof (struct AREANODE) *totalNode);

    int i = 0;
    for (i = 0; i < totalNode; i++) {
        areaTree[i].subArea = NULL;
        areaTree[i].nodeNum = 0;
    }
    return areaTree;
}

/********************
 *功能：输出分支信息到文件
 *输入：
 *输出：
 */
void showArea(char *fName, struct AREANODE *area) {

    printf("show area tree\n");
    FILE *fp;
    if ((fp = fopen(fName, "w")) == NULL) {
        printf("open file error!!\n");
        return;
    }

    int i = 0;
    int j = 0;
    for (i = 0; i < 360 * 180; i++) {

        if (area[i].nodeNum > 0) {
            j++;
            fprintf(fp, "%8d%5d%5d%8d", i + 1, i % 360, i / 360, area[i].nodeNum);
            struct SAMPLE *tmp = area[i].subArea;
            while (tmp) {
                fprintf(fp, "%15.8f", tmp->alpha);
                tmp = tmp->next;
            }
            fprintf(fp, "\n");
        }
    }
    printf("total number of area is:%d\n", j);
    fclose(fp);
}

/********************
 *功能：获得点所属分区
 *输入：点
 *输出：分支索引
 */
long getPointBranch(struct SAMPLE *point) {

    float alpha = (point->alpha - raMini);
    float delta = (point->delta - decMini);
    int x = (int) (alpha * factor);
    int y = (int) (delta * factor);

    return y * raNode + x;
}

//dec:angle, errorRadius: angle

float getAngleFromGreatCircle(double dec, double errorRadius) {
    double rst = acos((cos(errorRadius * ANG_TO_RAD) - pow(sin(dec * ANG_TO_RAD), 2)) / pow(cos(dec * ANG_TO_RAD), 2));
    return rst*RAD_TO_ANG;
}

/********************
 *功能：获得点point的可能搜索分支
 *输入：点
 *输出：搜索分支个数number，分支索引数组branch
 */
long * getPointSearchBranch(struct SAMPLE *point, long *number) {

    int height, width;

    float alpha = point->alpha;
    float delta = point->delta;

    float up = delta + searchRadius; //on north, up > down
    float down = delta - searchRadius; //on south, down > up
    float maxDec = 0.0;
    if (up > 0.0 && down > 0.0) {
        maxDec = up;
    } else if (up < 0.0 && down < 0.0) {
        maxDec = fabs(down);
    } else {
        if (fabs(up) > fabs(down))
            maxDec = fabs(up);
        else
            maxDec = fabs(down);
    }
    /*
        float raRadius = getAngleFromGreatCircle(maxDec, searchRadius);
     */
    int tIndex = ceil((maxDec - absDecMin) / searchRadius);
    float raRadius = raRadiusIndex[tIndex];

    float left = alpha - raRadius;
    float right = alpha + raRadius;

    //-zoneLength/3600
    if (up > 90.0) {
        up = 90.0;
        left = raMinf;
        right = raMaxf;
    } else if (down < -90.0) {
        down = -90.0;
        left = raMinf;
        right = raMaxf;
    }

    int indexUp = (up - decMini) * factor;
    int indexDown = (down - decMini) * factor;
    int indexLeft = (left - raMini) * factor;
    int indexRight = (right - raMini) * factor;

    if (indexUp >= decNode) indexUp = decNode - 1;
    if (indexDown < 0) indexDown = 0;
    if (indexRight >= raNode) indexRight = raNode - 1;
    if (indexLeft < 0) indexLeft = 0;

    height = abs(indexUp - indexDown) + 1;
    width = indexRight - indexLeft + 1;
    *number = height*width;
    long *branch = (long *) malloc(*number * sizeof (long));

    int i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            branch[i * width + j] = (indexDown + i) * raNode + indexLeft + j;
        }
    }

    return branch;
}

/********************
 *功能：将点添加point到分支branch
 *输入：
 *输出：
 */
void addPointToBranchSort(struct SAMPLE *point, struct AREANODE *branch) {

    point->next = NULL;
    struct SAMPLE *tmp = branch->subArea;
    if (tmp == NULL) {
        branch->subArea = point;
    } else {
        /**************************************/
        /*该部分是否考虑用二叉树或四叉树等其他优化查找的算法来实现*/
        /*当前按alpha的值从小到大排列*/
        /**/
        if (point->alpha < tmp->alpha) {
            branch->subArea = point;
            point->next = tmp;
        } else {
            struct SAMPLE *before = tmp;
            tmp = before->next;
            while ((tmp) && (point->alpha >= tmp->alpha)) { //当tmp的next为空时，tmp的下一个就是point的位置
                before = tmp;
                tmp = before->next;
            }
            before->next = point;
            point->next = tmp;
        }
    }
    branch->nodeNum = branch->nodeNum + 1;
}

void addPointToBranchNotSort(struct SAMPLE *point, struct AREANODE *branch) {

    point->next = NULL;
    struct SAMPLE *tmp = branch->subArea;
    if (tmp == NULL) {
        branch->subArea = point;
    } else {
        struct SAMPLE *before = tmp;
        tmp = before->next;
        before->next = point;
        point->next = tmp;
    }
    branch->nodeNum = branch->nodeNum + 1;
}

/********************
 *功能：将数据链表listHead添加到区域索引areaTree
 *输入：
 *输出：添加的总节点个数
 */
long addDataToTree(struct SAMPLE *head, struct AREANODE *areaTree) {

    long start, end;
    start = clock();

    struct SAMPLE *listHead = head;
    struct SAMPLE *tmp = listHead->next; //一般第一个节点为头结点，不存储数据，第二个节点才是第一个数据点
    long branch = 0;
    long i = 0;
    while (tmp) {
        listHead->next = tmp->next; //把tmp点从数据表中移除
        branch = getPointBranch(tmp); //获得tmp所属的树枝位置
        addPointToBranchNotSort(tmp, areaTree + branch); //把tmp点加入到树干的对应树枝中
        tmp = listHead->next; //取下一个点
        i++;
    }
    end = clock();
    if (showProcessInfo) {
        printf("totle point in index: %d\n", i);
        printf("time of init index is: %fs\n", (end - start)*1.0 / ONESECOND);
    }
    return i;
}

/********************
 *功能：返回总数据区域areaTree中有数据的区域的个数
 *输入：
 *输出：总数据区域中有数据的区域的个数
 */
int getTreeNodeNum(struct AREANODE *areaTree) {

    long num = 0;
    int i = 0;
    for (i = 0; i < 360 * 180; i++) {
        if (areaTree[i].nodeNum > 0) {
            num++;
        }
    }
    return num;
}

/********************
 *功能：将数值从角度转换为弧度
 *输入：角度
 *输出：弧度
 */
double angToRad(double angle) {
    return angle * ANG_TO_RAD;
}

/********************
 *功能：计算两个点的大圆距离
 *输入：两点
 *输出：距离
 */
double getGreatCircleDistance(struct SAMPLE *p1, struct SAMPLE *p2) {
    double rst = RAD_TO_ANG * acos(sin(ANG_TO_RAD * (p1->delta)) * sin(ANG_TO_RAD * (p2->delta)) +
            cos(ANG_TO_RAD * (p1->delta)) * cos(ANG_TO_RAD * (p2->delta)) * cos(ANG_TO_RAD * (fabs(p1->alpha - p2->alpha))));
    return rst;
}

/********************
 *功能：计算两个点得直线距离,p1 from raference catalog , p2 from sample catalog
 *输入：两点
 *输出：距离
 */
double getLineDistance(struct SAMPLE *p1, struct SAMPLE *p2) {
    double tmp = (pow(p1->pixx - p2->pixx1, 2) + pow(p1->pixy - p2->pixy1, 2));
    return sqrt(tmp);
}

/********************
 *功能：判断两个点是否匹配
 *输入：两点
 *输出：1匹配，0不匹配
 */
int isSimilarPoint(struct SAMPLE *p1, struct SAMPLE *p2) {
    double distance = 180.0 * acos(sin(angToRad(p1->delta)) * sin(angToRad(p2->delta)) +
            cos(angToRad(p1->delta)) * cos(angToRad(p2->delta)) * cos(angToRad(fabs(p1->alpha - p2->alpha)))) / PI;
    return distance < areaBox ? 1 : 0;
}

/********************
 *功能：在SAMPLE链表branch中寻找与point最匹配的点
 *输入：SAMPLE链表branch， POINTNOD点point
 *输出：误差error，目标掉goalPoint
 */
double searchSimilarPoint(struct SAMPLE *branch, struct SAMPLE *point, struct SAMPLE **goalPoint) {

    double error = areaBox;
    //double error = 1.0;

    /*branch按照alpha进行排序，可以匹配的区域（branch + - 0.016667）只是一个范围，
    进入这个范围之后，start=1，当出了这个范围时，说明后面的点已经不满足要求，直接break*/
    int start = 0;
    struct SAMPLE *tSample = branch;
    while (tSample) {

        //if((branch->alpha > minerror) && (branch->alpha < maxerror) ){
        /*
                if ((tSample->alpha + SUBAREA > point->alpha) && (tSample->alpha - SUBAREA < point->alpha)) {
         */
        start = 1;
        float distance = getGreatCircleDistance(tSample, point);
        //if(point->id == 9841 && tSample->id == 12569)
        //printf("diatance = %f\n", distance);
        //float distance = getLineDistance(branch, point);
        /**************************************/
        /*如果数据量大，上面建表时将采用二叉树等算法，这里采用对应的查找算法*/
        if (distance < error) {
            *goalPoint = tSample;
            error = distance;
            //break;		//现在是找到小于ERROR_GREAT_CIRCLE的点就结束循环，可以继续循环下去，找distance最小的点，有没有必要？
        }
        /*
                } else if (start == 1) { //当start==1时，说明point点已经从进入区域（branch + - 0.016667）后，又出来了
                    break;
                }
         */
        tSample = tSample->next;
    }
    return error;
}

/********************
 *功能：将sample dataB与reference areaTree表进行匹配
 *输入：sample dataB，reference areaTree
 *输出：匹配链表MATCHLIST
 */
void matchPoints(struct AREANODE *areaTree, struct SAMPLE *dataB) {

    long start, end;
    start = clock();
    /**/
    double error = areaBox;
    double minError = areaBox;
    //double error = 1.0;
    //double minError = 1.0;

    struct SAMPLE *samplePoint = NULL;
    struct SAMPLE *tmpPoint = NULL;
    struct SAMPLE *minPoint = NULL;
    samplePoint = dataB->next;
    long *branchIndex = NULL;

    int i = 0;
    int j = 0;
    int outData = 0;
    while (samplePoint) {
        if (samplePoint->alpha > raMaxi || samplePoint->alpha < raMini || samplePoint->delta > decMaxi || samplePoint->delta < decMini) {
            samplePoint->reference = NULL;
            samplePoint->crossid = -1;
            samplePoint->error = 100;
            samplePoint = samplePoint->next;
            outData++;
            continue;
        }

        long numArea = 0;
        branchIndex = getPointSearchBranch(samplePoint, &numArea);

        //error = minError = 1.0;
        error = minError = areaBox;

        minPoint = NULL;
        tmpPoint = NULL;

        for (i = 0; i < numArea; i++) {
            struct SAMPLE *branch = areaTree[branchIndex[i]].subArea;
            error = searchSimilarPoint(branch, samplePoint, &tmpPoint);
            if (minError > error && tmpPoint != NULL) {
                minError = error;
                minPoint = tmpPoint;
            }
        }
        if (minPoint) {
            samplePoint->reference = minPoint;
            samplePoint->crossid = minPoint->id;
            samplePoint->error = minError;
        } else {
            samplePoint->reference = NULL;
            samplePoint->crossid = -1;
            samplePoint->error = 100;
        }
        samplePoint = samplePoint->next;
        j++;
    }

    end = clock();
    if (showProcessInfo) {
        printf("time of cross match is: %fs\n", (end - start)*1.0 / ONESECOND);
        printf("out data: %d\n", outData);
    }

}


void freeList(void *link) {
    long next;
    if (link == NULL) return;
    while (link != NULL) {
        next = *((long*) link);
        free(link);
        link = (void*) next;
    }
    return;
}


void writeNotMatchedToFile(char *fileName, struct SAMPLE * sample) {

    if (showProcessInfo) {
        printf("start write file!\n");
    }
    long start, end;
    start = clock();

    if (sample == NULL) return;
    FILE *fp = NULL;

    char *nameBuf = (char*) malloc(256);
    char *tmpBuf = (char*) malloc(256);

    char *pos1, *pos2;
    pos1 = strrchr(fileName, '\\');
    if (pos1 == NULL) {
        pos1 = strrchr(fileName, '/');
        if (pos1 == NULL) {
            pos1 = fileName;
        } else {
            pos1 = pos1 + 1;
        }
    } else {
        pos1 = pos1 + 1;
    }
    pos2 = strrchr(fileName, '.');
    if (pos2 == NULL)
        pos2 = fileName + strlen(fileName);
    strncpy(nameBuf, pos1, pos2 - pos1);
    nameBuf[pos2 - pos1] = '\0';
    /*
        printf("\nnew file name= %s\n", nameBuf);
     */
    sprintf(tmpBuf, "log_%s.txt", nameBuf);
    /*
        printf("new file name= %s\n\n", tmpBuf);
     */


    if ((fp = fopen(tmpBuf, "w")) == NULL) {
        printf("open file %s error!\n", tmpBuf);
        return;
    }
    fprintf(fp, "%10s%15s%15s%10s%15s%15s%15s", "ida", "alpha", "delta", "idb", "alpha", "delta", "error\n");

    struct SAMPLE *tSample = sample->next;
    long i = 0;
    while (tSample) {
        if (tSample->reference == NULL) {
            fprintf(fp, "%10d %15.7f %15.7f %10d %15.7f %15.7f %15.7f\n",
                    tSample->id, tSample->alpha, tSample->delta,
                    -1, -1.0, -1.0, -1.0);
            i++;
        } else if (tSample->error >= areaBox) {
            fprintf(fp, "%10d %15.7f %15.7f %10d %15.7f %15.7f %15.7f\n",
                    tSample->id, tSample->alpha, tSample->delta,
                    tSample->reference->id, tSample->reference->alpha, tSample->reference->delta,
                    tSample->error);
            i++;
        }
        /*
                if(tSample->error<areaBox){
                    fprintf(fp, "%10d %15.7f %15.7f %10d %15.7f %15.7f %15.7f\n", 
                        tSample->id, tSample->alpha, tSample->delta, 
                        tSample->reference->id, tSample->reference->alpha, tSample->reference->delta, 
                        tSample->error);
                    i++;
                }
         */
        tSample = tSample->next;
    }
    fclose(fp);

    end = clock();
    if (showProcessInfo) {
        printf("total unmatched stars:%d\n", i);
        printf("the unmatched stars was written to file %s use time: %fs\n", tmpBuf, (end - start)*1.0 / ONESECOND);
    }
    free(tmpBuf);
    free(nameBuf);
}

void writeToTXT(char *fileName, struct SAMPLE * sample) {

    if (showProcessInfo) {
        printf("start write file!\n");
    }

    if (sample == NULL) return;
    FILE *fp = NULL;

    char *nameBuf = (char*) malloc(256);
    char *tmpBuf = (char*) malloc(256);

    char *pos1, *pos2;
    pos1 = strrchr(fileName, '\\');
    if (pos1 == NULL) {
        pos1 = strrchr(fileName, '/');
        if (pos1 == NULL) {
            pos1 = fileName;
        } else {
            pos1 = pos1 + 1;
        }
    } else {
        pos1 = pos1 + 1;
    }
    pos2 = strrchr(fileName, '.');
    if (pos2 == NULL)
        pos2 = fileName + strlen(fileName);
    strncpy(nameBuf, pos1, pos2 - pos1);
    nameBuf[pos2 - pos1] = '\0';
    /*
        printf("\nnew file name= %s\n", nameBuf);
     */
    sprintf(tmpBuf, "%s.txt", nameBuf);
    /*
        printf("new file name= %s\n\n", tmpBuf);
     */


    if ((fp = fopen(tmpBuf, "w")) == NULL) {
        printf("open file %s error!\n", tmpBuf);
        return;
    }
    fprintf(fp, "#%10s%15s%15s", "ida", "alpha", "delta\n");

    struct SAMPLE *tSample = sample->next;
    long i = 0;
    while (tSample) {
        fprintf(fp, "%10d %15.7f %15.7f\n", tSample->id, tSample->alpha, tSample->delta);
        i++;
        tSample = tSample->next;
    }
    fclose(fp);

    if (showProcessInfo) {
        printf("write file %s total stars:%d\n", tmpBuf, i);
    }
    free(tmpBuf);
    free(nameBuf);
}

struct SAMPLE * matchAll(struct SAMPLE *dataA, struct SAMPLE * dataB) {

    struct SAMPLE *ta = NULL;
    struct SAMPLE *tb = dataB->next;

    struct SAMPLE *tb2, *dataB2Head;
    struct SAMPLE *dataB2 = (struct SAMPLE*) malloc(sizeof (struct SAMPLE));
    memcpy(dataB2, dataB, sizeof (struct SAMPLE));
    dataB2->next = NULL;
    dataB2Head = dataB2;

    long i = 0;
    while (tb) {
        ta = dataA->next;
        tb2 = (struct SAMPLE*) malloc(sizeof (struct SAMPLE));
        memcpy(tb2, tb, sizeof (struct SAMPLE));
        tb2->error = 100.0;
        tb2->next = NULL;
        dataB2->next = tb2;
        dataB2 = tb2;
        while (ta) {
            float distance = getGreatCircleDistance(ta, tb2);
            if (tb2->error > distance) {
                tb2->error = distance;
                tb2->reference = ta;
            }
            ta = ta->next;
            i++;
        }
        tb = tb->next;
    }
    return dataB2Head;
}

void writeDifferentToFile(char *fileName, struct SAMPLE *crossRst, struct SAMPLE * zoneRst) {

    if (showProcessInfo) {
        printf("start write different to file!\n");
    }
    long start, end;
    start = clock();

    if (crossRst == NULL || zoneRst == NULL) return;
    FILE *fp = NULL;

    char *nameBuf = (char*) malloc(256);
    char *tmpBuf = (char*) malloc(256);


    char *pos1, *pos2;
    pos1 = strrchr(fileName, '\\');
    if (pos1 == NULL) {
        pos1 = strrchr(fileName, '/');
        if (pos1 == NULL) {
            pos1 = fileName;
        } else {
            pos1 = pos1 + 1;
        }
    } else {
        pos1 = pos1 + 1;
    }
    pos2 = strrchr(fileName, '.');
    if (pos2 == NULL)
        pos2 = fileName + strlen(fileName);
    strncpy(nameBuf, pos1, pos2 - pos1);
    nameBuf[pos2 - pos1] = '\0';
    /*
        printf("\nnew file name= %s\n", nameBuf);
     */
    sprintf(tmpBuf, "cross_different_zone_%s.txt", nameBuf);
    /*
        printf("new file name= %s\n\n", tmpBuf);
     */


    if ((fp = fopen(tmpBuf, "w")) == NULL) {
        printf("open file %s error!\n", tmpBuf);
        return;
    }
    fprintf(fp, "#%10s%15s%15s%10s%15s%15s%15s", "ida", "alpha", "delta", "idb", "alpha", "delta", "error\n");

    struct SAMPLE *tCrossRst = crossRst->next;
    struct SAMPLE *tZoneRst = zoneRst->next;
    long crossMatched = 0;
    long zoneMatched = 0;
    while (tCrossRst) {
        if (tZoneRst->id == tCrossRst->id) {
            if (tCrossRst->error < areaBox) {
                if (tZoneRst->error >= areaBox) { //find the cross method matched but the zone method not matched, and output to file
                    fprintf(fp, "%10d %15.7f %15.7f %10d %15.7f %15.7f %15.7f\n",
                            tCrossRst->id, tCrossRst->alpha, tCrossRst->delta,
                            tCrossRst->reference->id, tCrossRst->reference->alpha, tCrossRst->reference->delta,
                            tCrossRst->error);
                    zoneMatched++;
                }
                crossMatched++;
            }
        } else {
            printf("copy sample data error, the sample data1(cross result) is different sampel data2(zone resutlt)! please check!\n");
        }
        tCrossRst = tCrossRst->next;
        tZoneRst = tZoneRst->next;
    }
    fclose(fp);

    end = clock();
    printf("cross method matched starts: %d\n", crossMatched);
    printf("zone method omitted stars: %d \n", zoneMatched);
    printf("the omitted star was written to file %s use time: %fs\n", tmpBuf, (end - start)*1.0 / ONESECOND);

    free(tmpBuf);
    free(nameBuf);
}

void showPartitionInfo(char *fName, struct AREANODE * area) {

    printf("show area tree\n");
    FILE *fp;
    if ((fp = fopen(fName, "w")) == NULL) {
        printf("open file error!!\n");
        return;
    }

    int i = 0;
    int j = 0;
    for (i = 0; i < planeZoneX * planeZoneY; i++) {

        if (area[i].nodeNum > 0) {
            j++;
            fprintf(fp, "%8d%5d%5d%8d", i + 1, i % planeZoneX, i / planeZoneX, area[i].nodeNum);
            struct SAMPLE *tmp = area[i].subArea;
            while (tmp) {
                fprintf(fp, "%15.8f%15.8f", tmp->pixx, tmp->pixy);
                tmp = tmp->next;
            }
            fprintf(fp, "\n");
        }
    }
    printf("total number of area is:%d\n", j);
    fclose(fp);
}

double searchSimilarPointPlane(struct SAMPLE *branch, struct SAMPLE *point, struct SAMPLE **goalPoint) {

    double error = areaBox;

    struct SAMPLE *tSample = branch;
    while (tSample) {
        float distance = getLineDistance(tSample, point);
        if (distance < error) {
            *goalPoint = tSample;
            error = distance;
        }
        tSample = tSample->next;
    }
    return error;
}

struct SAMPLE *readStarFile(char *fName, int &starNum) {

  FILE *fp = fopen(fName, "r");
  if (fp == NULL) {
    return NULL;
  }

  starNum = 0;
  float alpha=0, delta=0;
  struct SAMPLE *starList = (struct SAMPLE *) malloc(sizeof (struct SAMPLE));
  starList->next = NULL;
  struct SAMPLE *tStar = NULL;
  struct SAMPLE *nextStar = NULL;
  char line[MaxStringLength];

  while (fgets(line, MaxStringLength, fp) != NULL) {
    if (2 == sscanf(line, "%f%f%*s", &alpha, &delta)) {
      nextStar = (struct SAMPLE *) malloc(sizeof (struct SAMPLE));
      nextStar->id = starNum;
      nextStar->alpha = alpha;
      nextStar->delta = delta;
      nextStar->next = NULL;
      if (NULL == starList->next) {
        starList->next = nextStar;
        tStar = nextStar;
      } else {
        tStar->next = nextStar;
        tStar = nextStar;
      }
      starNum++;
    }
  }

  printf("%s\t%d stars\n", fName, starNum);
  
  return starList;
}

void writeToFile(char *fileName, struct SAMPLE * sample) {

    if (showProcessInfo) {
        printf("start write file!\n");
    }

    if (sample == NULL) return;
    FILE *fp = NULL;

    if ((fp = fopen(fileName, "w")) == NULL) {
        printf("open file %s error!\n", fileName);
        return;
    }
    fprintf(fp, "#%10s%15s%15s%10s%15s%15s", "sid", "sra", "sdec", "rid", "rra", "rdec", "error\n");

    struct SAMPLE *tSample = sample->next;
    long i = 0;
    while (tSample) {
      if(tSample->reference != NULL){
        fprintf(fp, "%10d %15.7f %15.7f%10d %15.7f %15.7f\n", tSample->id, tSample->alpha, tSample->delta, tSample->reference->id, tSample->reference->alpha, tSample->reference->delta, tSample->error);
        i++;
      }
        tSample = tSample->next;
    }
    fclose(fp);

    if (1) {
        printf("write file %s total stars:%d\n", fileName, i);
    }
}


void mainSphere(char *inFile1, char *inFile2, char *outFile) {

    //printf("starting corss match...\n");

    long start, end;
    start = clock();
    struct SAMPLE *dataB2;
    
    int numA = 0, numB = 0;

    struct SAMPLE *dataA = readStarFile(inFile1, numA);
    //writeToTXT(inFile1, dataA);
    if (dataA == NULL) return;

    struct SAMPLE *dataB = readStarFile(inFile2, numB);
    //writeToTXT(inFile2, dataB);
    if (dataB == NULL) return;

    if (useCross == 1) {
        dataB2 = matchAll(dataA, dataB);
    }

    struct AREANODE *areaTree = initAreaNode(dataA);

    addDataToTree(dataA, areaTree);
    //showArea("data/showArea.txt", areaTree) ;
    matchPoints(areaTree, dataB);

    writeToFile(outFile, dataB);

    if (useCross == 1) {
        writeDifferentToFile(inFile2, dataB2, dataB);
        freeList(dataB2);
    }
    if (showProcessInfo) {
        writeNotMatchedToFile(inFile2, dataB);
    }
    end = clock();
    if (showProcessInfo) {
        printf("total time is: %fs\n", (end - start)*1.0 / ONESECOND);
        printAreaInfo();
    }
    printf("total time is: %fs\n", (end - start)*1.0 / ONESECOND);
    printf("corss match done!\n");

    free(areaTree);
    freeList(dataA);
    freeList(dataB);
    if (useCross == 1) {
        freeList(dataB2);
    }

}