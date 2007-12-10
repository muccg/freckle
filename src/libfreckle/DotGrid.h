#ifndef _DOTGRID_H_
#define _DOTGRID_H_

#include "DotStore.h"


#include <assert.h>

class DotGrid
{
private:
	int width;
	int height;

	int *data;			// where we store the grid data

	// for internal class use (just to make more robust and bounds test)
	inline void SetPoint(int x, int y, int point)
	{
		assert(data);
		assert(x>=0 && x<width);
		assert(y>=0 && y<height);
		data[y*width+x]=point;
	}

	inline int GetPoint(int x, int y) const
	{
		assert(data);
		assert(x>=0 && x<width);
		assert(y>=0 && y<height);
		return data[y*width+x];
	}

	inline int GetData(int pos) const
	{
		assert(data);
		assert(pos>=0 && pos<width*height);
		return data[pos];
	}

public:
	DotGrid();
	~DotGrid();

	void Create(int xsize, int ysize);
	void Destroy();

	// some helpers
	inline int GetWidth() const
	{
		return width;
	}

	inline int GetHeight() const
	{
		return height;
	}

	inline int GetSize() const
	{
		return width*height;
	}	

	//
	// \brief helper function to calculate maximum array value (for scaling)
	//
	// \param 
	inline int GetMax() const
	{
		assert(data);
		int max=-999999;
		for(int pos=0; pos<width*height; pos++)
			if(data[pos]>max)
				max=data[pos];
		return max;
	} 

	//
	// \brief helper function to calculate maximum array value (for scaling)
	//
	// \param 
	inline int GetMin() const
	{
		assert(data);
		int min=99999999;
		for(int pos=0; pos<width*height; pos++)
			if(data[pos]<min)
				min=data[pos];
		return min;
	} 

	// turn the averaged grid into a luminance image string for the higher level language
	unsigned char *ToString() const;

	// \brief Taking data from the specified source, calculate an aveage grid using the window (x1,y1) - (x2,y2) as the source area
	//
	// \param source The DotStore from which to take our dot data
	// \param x1 the x position of the top left corner of the source window
	// \param y1 the y position of the top left corner of the source window
	// \param x2 the x position of the bottom right corner of the source window
	// \param y2 the y position of the bottom right corner of the source window
	// \param scale how much to shrink the window by in averaging. A value of 10 means make the grid 1/10th of the original window size. scale must be an even divisor of both the width and the height of the souce window
	// \param window the window size used during the DotStore calculations: TODO: grab this value from the dotstore itself
	void CalculateGrid(DotStore *source, double x1, double y1, double x2, double y2, double scale, int window);

	// Add the values from another grid. Grid must have the same dimensions. Used to combine forward and reverse plots
	void AddInplace(DotGrid *second);

	// flip this grid upside down
	void FlipInplace();

};

#endif

