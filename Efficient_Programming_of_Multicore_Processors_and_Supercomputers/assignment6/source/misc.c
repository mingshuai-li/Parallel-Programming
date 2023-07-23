/*
 * misc.c
 *
 * Helper functions for
 * - initialization
 * - finalization,
 * - writing out a picture
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#include "heat.h"

/*
 * Initialize the iterative solver
 * - allocate memory for matrices
 * - set boundary conditions according to configuration
 */
int initialize( algoparam_t *param )
{
    int i, j;
    double dist;

    // total number of points (including border)
    const int np = param->act_res + 2;

    //
    // allocate memory
    //
    (param->u)     = (double*)malloc( sizeof(double)* (param->numberOfRowsIncludingGhost)*(param->numberOfColumnsIncludingGhost) );
    (param->uhelp) = (double*)malloc( sizeof(double)* (param->numberOfRowsIncludingGhost)*(param->numberOfColumnsIncludingGhost) );
    (param->uvis)  = (double*)calloc( sizeof(double),
				      (param->visres+2) *
				      (param->visres+2) );
					  
	#pragma omp parallel for schedule(static)
    for (i=0;i<param->numberOfRowsIncludingGhost;i++){
    	for (j=0;j<param->numberOfColumnsIncludingGhost;j++){
    		param->u[i*param->numberOfColumnsIncludingGhost+j]=0.0;
			param->uhelp[i*param->numberOfColumnsIncludingGhost+j]=0.0;
    	}
    }

    if( !(param->u) || !(param->uhelp) || !(param->uvis) )
    {
	fprintf(stderr, "Error: Cannot allocate memory\n");
	return 0;
    }

    for( i=0; i<param->numsrcs; i++ )
    {
	/* top row */
	if(param->startRowIndex==1){
		for( j=0; j<param->numberOfColumns; j++ )
		{
			dist = sqrt(pow((double)(param->startColumnIndex+j)/(double)(np-1) - param->heatsrcs[i].posx, 2)+pow(param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[j+1] +=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[j+1]=(param->u)[j+1];
			}		
		}
		// west-north corner
		// Although the corner points can be processed together in the previous loop, 
		// we do not want to redundantly store values that we do not need on certain points in the ghost area, 
		// even though it does not affect the computation results.
		if(param->startColumnIndex==1){
			dist = sqrt(pow(param->heatsrcs[i].posx, 2)+pow(param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[0] +=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[0]=(param->u)[0];
			}	
		}
		// east-north corner
		if(param->endColumnIndex==(np-2)){
			dist = sqrt(pow(1.0-param->heatsrcs[i].posx, 2)+pow(param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[param->numberOfColumns+1] +=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[param->numberOfColumns+1]=(param->u)[param->numberOfColumns+1];
			}			
		}
	}

	/* bottom row */
	if(param->endRowIndex==(np-2)){
		for( j=0; j<param->numberOfColumns; j++ )
		{
			dist = sqrt(pow((double)(param->startColumnIndex+j)/(double)(np-1) - param->heatsrcs[i].posx, 2)+pow(1.0-param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)+j+1] +=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)+j+1]=(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)+j+1];
			}		
		}
		// west-south corner
		if(param->startColumnIndex==1){
			dist = sqrt(pow(param->heatsrcs[i].posx, 2)+pow(1.0-param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)]+=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)]=(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+1)];
			}	
		}
		// east-south corner
		if(param->endColumnIndex==(np-2)){
			dist = sqrt(pow(1.0-param->heatsrcs[i].posx, 2)+pow(1.0-param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+2)-1] +=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+2)-1]=(param->u)[param->numberOfColumnsIncludingGhost*(param->numberOfRows+2)-1];
			}			
		}
	}	

	/* leftmost column */
	if(param->startColumnIndex==1){
		for( j=0; j<param->numberOfRows; j++ )
		{
			dist = sqrt( pow(param->heatsrcs[i].posx, 2)+
				pow((double)(param->startRowIndex+j)/(double)(np-1) -
					param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[ (j+1)*param->numberOfColumnsIncludingGhost ]+=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[ (j+1)*param->numberOfColumnsIncludingGhost ] = (param->u)[ (j+1)*param->numberOfColumnsIncludingGhost ];
			}
		}
	}

	/* rightmost column */
	if(param->endColumnIndex==(np-2)){
		for( j=0; j<param->numberOfRows; j++ )
		{
			dist = sqrt( pow(1.0-param->heatsrcs[i].posx, 2)+
				pow((double)(param->startRowIndex+j)/(double)(np-1) -
					param->heatsrcs[i].posy, 2));
			if( dist <= param->heatsrcs[i].range )
			{
			(param->u)[ (j+2)*param->numberOfColumnsIncludingGhost-1 ]+=
				(param->heatsrcs[i].range-dist) /
				param->heatsrcs[i].range *
				param->heatsrcs[i].temp;
			(param->uhelp)[ (j+2)*param->numberOfColumnsIncludingGhost-1 ] = (param->u)[ (j+2)*param->numberOfColumnsIncludingGhost-1 ];
			}
		}
	}

    }

    return 1;
}

/*
 * free used memory
 */
int finalize( algoparam_t *param )
{
    if( param->u ) {
	free(param->u);
	param->u = 0;
    }

    if( param->uhelp ) {
	free(param->uhelp);
	param->uhelp = 0;
    }

    if( param->uvis ) {
	free(param->uvis);
	param->uvis = 0;
    }

    return 1;
}


/*
 * write the given temperature u matrix to rgb values
 * and write the resulting image to file f
 */
void write_image( FILE * f, double *u,
		  unsigned sizex, unsigned sizey )
{
    // RGB table
    unsigned char r[1024], g[1024], b[1024];
    int i, j, k;

    double min, max;

    j=1023;

    // prepare RGB table
    for( i=0; i<256; i++ )
    {
	r[j]=255; g[j]=i; b[j]=0;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=255-i; g[j]=255; b[j]=0;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=0; g[j]=255; b[j]=i;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=0; g[j]=255-i; b[j]=255;
	j--;
    }

    min=DBL_MAX;
    max=-DBL_MAX;

    // find minimum and maximum
    for( i=0; i<sizey; i++ )
    {
	for( j=0; j<sizex; j++ )
	{
	    if( u[i*sizex+j]>max )
		max=u[i*sizex+j];
	    if( u[i*sizex+j]<min )
		min=u[i*sizex+j];
	}
    }


    fprintf(f, "P3\n");
    fprintf(f, "%u %u\n", sizex, sizey);
    fprintf(f, "%u\n", 255);

    for( i=0; i<sizey; i++ )
    {
	for( j=0; j<sizex; j++ )
	{
	    k=(int)(1024.0*(u[i*sizex+j]-min)/(max-min));
	    fprintf(f, "%d %d %d  ", r[k], g[k], b[k]);
	}
	fprintf(f, "\n");
    }
}

int coarsen( double *uold, unsigned oldx, unsigned oldy ,
		int startColumnIndex, int startRowIndex, int columnsGhost, int rowsGhost,
	    double *unew, unsigned newx, unsigned newy )
{
    int i, j, k, l, ii, jj;

	int endColumnIndex = startColumnIndex + columnsGhost - 3;
	int endRowIndex = startRowIndex + rowsGhost - 3;
    int stopx = newx;
    int stopy = newy;
    float temp;
    float stepx = (float)oldx/(float)newx;
    float stepy = (float)oldy/(float)newy;

    if (oldx<newx){
	 stopx=oldx;
	 stepx=1.0;
    }
    if (oldy<newy){
     stopy=oldy;
     stepy=1.0;
    }

    //printf("oldx=%d, newx=%d\n",oldx,newx);
    //printf("oldy=%d, newy=%d\n",oldy,newy);
    //printf("rx=%f, ry=%f\n",stepx,stepy);
    // NOTE: this only takes the top-left corner,
    // and doesnt' do any real coarsening
	
	// We handle each section separately to avoid accumulating redundant values.
	/* contribution from inner part without boundary*/
    for( i=0; i<stopy; i++ ){
       ii=stepy*i;
       for( j=0; j<stopx; j++ ){
          jj=stepx*j;
          for ( k=0; k<stepy; k++ ){
	       	for ( l=0; l<stepx; l++ ){
	       		if (ii+k>=startRowIndex && jj+l>=startColumnIndex 
				&& ii+k<=endRowIndex && jj+l<=endColumnIndex){
		           unew[i*newx+j] += uold[(ii+k+1-startRowIndex)*columnsGhost+(jj+l+1-startColumnIndex)] ;
				}
			}
	      }    
       }
    }

	/* north boundary without corner */
	if (startRowIndex==1){
		for( j=0; j<stopx; j++ ){
			jj=stepx*j;
			for ( l=0; l<stepx; l++ ){
				if (jj+l>=startColumnIndex && jj+l<=endColumnIndex){
					unew[j] += uold[jj+l+1-startColumnIndex] ;
				}
			} 
		}
	}

	/* south boundary without corner*/
	if (endRowIndex==(oldy-2)){
		for( j=0; j<stopx; j++ ){
			jj=stepx*j;
			for ( l=0; l<stepx; l++ ){
				if (jj+l>=startColumnIndex && jj+l<=endColumnIndex){
					unew[(stopy-1)*newx+j] += uold[(oldy-startRowIndex)*columnsGhost+(jj+l+1-startColumnIndex)] ;
				}
			} 
		}
	}

	/* west boundary without corner */
	if(startColumnIndex==1){
		for( i=0; i<stopy; i++ ){
		ii=stepy*i;
			for ( k=0; k<stepy; k++ ){
				if (ii+k>=startRowIndex && ii+k<=endRowIndex){
					unew[i*newx] += uold[(ii+k+1-startRowIndex)*columnsGhost] ;
				}
			}
		}
	}
	/* east boundary without corner */
	if(endColumnIndex==(oldx-2)){
		for( i=0; i<stopy; i++ ){
		ii=stepy*i;
			for ( k=0; k<stepy; k++ ){
				if (ii+k>=startRowIndex && ii+k<=endRowIndex){
					unew[i*newx+(stopx-1)] += uold[(ii+k+1-startRowIndex)*columnsGhost+(oldx-startColumnIndex)] ;
				}
			}
		}	
	}
	
	/* west-north corner */
	if(startRowIndex==1 && startColumnIndex==1){
		unew[0] += uold[0];
	}
   
	/* west-south corner */
	if(endRowIndex==(oldy-2) && startColumnIndex==1){
		unew[(stopy-1)*newx] += uold[(rowsGhost-1)*columnsGhost];
	}

	/* east-north corner */
	if(startRowIndex==1 && endColumnIndex==(oldx-2)){
		unew[stopx-1] += uold[columnsGhost-1];
	}	    

	/* east-south corner */
	if(endRowIndex==(oldy-2) && endColumnIndex==(oldx-2)){
		unew[(stopy-1)*newx+(stopx-1)] += uold[rowsGhost*columnsGhost-1];
	}

	/* average */
	float tmp = 1/(stepx*stepy);
	for (i = 0; i < newy; ++i){
		for (j = 0; j < newx; ++j){
			unew[i*newx+j] *= tmp;
		}
	}

  return 1;
}
