#include "DotGrid.h"
#include <memory.h>
#include <assert.h>

DotGrid::DotGrid()
{
	width=height=0;
	data=NULL;
}

DotGrid::~DotGrid()
{
	Destroy();
}

void DotGrid::Create(int xsize, int ysize)
{
	data=new int [xsize*ysize];

	memset(data,0,xsize*ysize*sizeof(int));

	assert(data);

	width=xsize;
	height=ysize;
}

void DotGrid::Destroy()
{
	if(data)
		delete data;

	data=NULL;
}

unsigned char *DotGrid::ToString() const
{
	unsigned char *out=new unsigned char [width*height];
	assert(out);	

	// normalise the histogram
	int numpixels=width*height;
	int *histogram=CalculateHistogram();

	assert(histogram);

	// print our histogram
// 	for(int i=0; i<GetMax()+1;i++)
// 		printf("%d => %d\n",i,histogram[i]);

	for(int pos=0; pos<width*height; pos++)
		out[pos]=(unsigned char)(255.-255.*(double)numpixels*((double)histogram[data[pos]]-(double)histogram[0])/(((double)numpixels-(double)histogram[0])*(double)numpixels));

	// free histogram
	delete histogram;

	return out;
}

int *DotGrid::CalculateHistogram() const
{
	int max=GetMax();

	int *histogram= new int[max+1];
	assert(histogram);
	memset(histogram, 0, sizeof(int)*(max+1));

	for(int pos=0; pos<width*height; pos++)
		histogram[data[pos]]++;

	//cumulative histogram
	for(int i=1; i<max+1;i++)
		histogram[i]+=histogram[i-1];

	return histogram;
}

// Calculate a sum from a window on the source dotstore
void DotGrid::CalculateGrid(DotStore *source, double x1, double y1, double x2, double y2, double scale, int window)
{
	double xsize=x2-x1;
	double ysize=y2-y1;

	// the start and end values for the window should be in the proper order. top left to bottom right
	assert(xsize>0);
	assert(ysize>0);

	// check we divide evenly
// 	assert( (xsize / scale) - (double)(int)(xsize / scale) == 0.0 );
// 	assert( (ysize / scale) - (double)(int)(ysize / scale) == 0.0 );

	// the size of the averaging window and the remaining strips along the side and bottom
	int numx=(int)(xsize/scale);
	int numy=(int)(ysize/scale);

	// make sure we haven't already got a grid
	assert(!data);

	Create(numx,numy);

	for(int y=0; y<numy; y++)
		for(int x=0; x<numx; x++)
			SetPoint(x,y,source->CountAreaMatches(x*scale+x1, y*scale+y1, (x+1)*scale+x1, (y+1)*scale+y1, window) );
}

// Add the values from another grid. Grid must have the same dimensions. Used to combine forward and reverse plots
void DotGrid::AddInplace(DotGrid *second)
{
	assert(second->GetWidth()==GetWidth());
	assert(second->GetHeight()==GetHeight());

	// robust version
// 	for(int y=0; y<GetHeight(); y++)
// 		for(int x=0; x<GetWidth(); x++)
// 			SetPoint(x,y,GetPoint(x,y)+second->GetPoint(x,y));

	// fast version
	for(int pos=0; pos<width*height; pos++)
		data[pos]+=second->GetData(pos);
}

// flip this grid upside down
void DotGrid::FlipInplace()
{
	int *store=new int[width];
	for(int y=0; y<height/2; y++)
	{
		//copy row into store
		memcpy(store, &data[y*width], width*sizeof(int));
		
		//copy end row into head row
		memcpy(&data[y*width], &data[(height-1-y)*width], width*sizeof(int));

		// copy store into end row
		memcpy(&data[(height-1-y)*width], store, width*sizeof(int));
	}
	delete store;
}


