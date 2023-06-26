/* commonIntervals2.c */
/* (spikes, windowSize, numCommon, numCells)*/
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *spikes, *indexing, *out, *tempOutput, *cellOccurrances, windowSize, currMin;
    int i, j, k, lengthSpikes, row, numCommon, numCells, counter;
    
    /* Ensure inputs are correctly attributed*/
    if ( mxGetN(prhs[0]) != 2 ) mexErrMsgTxt("spikes must be an M x 2 matrix. The first column denotes interval start times, second column denotes interval end times.");
    
    /* Allocate values*/
    spikes = mxGetPr(prhs[0]);
    windowSize = (double)*mxGetPr(prhs[1]);
    numCommon = (int)*mxGetPr(prhs[2]);
    numCells = (int)*mxGetPr(prhs[3]);
    lengthSpikes = mxGetM(prhs[0]); /* Number of rows*/
    tempOutput = mxGetPr(mxCreateNumericMatrix(lengthSpikes, 2, mxDOUBLE_CLASS, mxREAL));
    cellOccurrances = mxGetPr(mxCreateNumericMatrix(numCells, 2, mxDOUBLE_CLASS, mxREAL));
    
    row = 0;
    j = 0;
    i = 0;
    /* The algorithm*/
    while (i < lengthSpikes-1)
    {
        j++;
        cellOccurrances = mxGetPr(mxCreateNumericMatrix(numCells, 1, mxDOUBLE_CLASS, mxREAL));
        cellOccurrances[((int) spikes[i+lengthSpikes])-1] = spikes[i];

        while (((spikes[j] - spikes[i]) < windowSize) && (j < lengthSpikes))
        {

            cellOccurrances[((int) spikes[j+lengthSpikes])-1] = spikes[j];
            if (spikes[j+lengthSpikes] == spikes[i+lengthSpikes])
            {
                cellOccurrances[((int) spikes[i+lengthSpikes])-1] = spikes[j];
            }
            
            counter = 0;
            for (k = 0; k < numCells; k++)
            {
                if (cellOccurrances[k] > 0)
                {
                    counter++;
                }
            }
            if (counter >= numCommon)
            {
                currMin = cellOccurrances[0];
                for (k=0; k<numCells;k++)
                {
                    if (cellOccurrances[k] > 0)
                    {
                        if (cellOccurrances[k]<currMin)
                        {
                            currMin = cellOccurrances[k];
                        }
                    }
                }
                tempOutput[row] = currMin;
                tempOutput[row+lengthSpikes] = spikes[j];
                row++;
                i = j-1;
                break;
            }
            j++;            
        }
        i++;
        j=i;
    }
    if (row != 0)
    {
        plhs[0] = mxCreateNumericMatrix(row, 2, mxDOUBLE_CLASS, mxREAL);
        out = mxGetPr(plhs[0]);
        
        for (i = 0; i < row; ++i)
        {
            out[i] = tempOutput[i];
            out[i+row] = tempOutput[i+lengthSpikes];
        }
    }
    else
    {
        plhs[0] = mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);
    }
}