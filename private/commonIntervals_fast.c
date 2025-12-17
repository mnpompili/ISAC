/* commonIntervals_fast.c */
/* (spikes, windowSize, numCells, minCellID)*/
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *spikes, *out, *tempOutput, *cellOccurrances, windowSize, currMin;
    int i, j, k, lengthSpikes, row, numCells, minCellID, counter, exitLoop, oldEnd;
    
    /* Ensure inputs are correctly attributed*/
    if (mxGetN(prhs[0]) != 2) mexErrMsgTxt("spikes must be an M x 2 matrix. The first column denotes interval start times, second column denotes interval end times.");
    
    /* Allocate values*/
    spikes = mxGetPr(prhs[0]);
    windowSize = (double)*mxGetPr(prhs[1]);
    numCells = (int)*mxGetPr(prhs[2]);
    minCellID = (int)*mxGetPr(prhs[3]);
    lengthSpikes = mxGetM(prhs[0]); /* Number of rows*/
    tempOutput = mxGetPr(mxCreateNumericMatrix(lengthSpikes, 2, mxDOUBLE_CLASS, mxREAL));
    cellOccurrances = mxGetPr(mxCreateNumericMatrix(numCells, 2, mxDOUBLE_CLASS, mxREAL));

    row = 0;
    j = 0;
    i = 0;
    oldEnd = 0;
    /* The algorithm*/
    for (i = 0; i<lengthSpikes; i++)
    {   
        if (spikes[i+lengthSpikes] == minCellID)
        {
            cellOccurrances = mxGetPr(mxCreateNumericMatrix(numCells, 1, mxDOUBLE_CLASS, mxREAL));
            cellOccurrances[minCellID-1] = spikes[i];
            j = i-1;
            counter = 1;
            while (((spikes[i] - spikes[j]) < windowSize) && (j > oldEnd))             
            {
                if (cellOccurrances[((int) spikes[j+lengthSpikes])-1] == 0)
                {
                    cellOccurrances[((int) spikes[j+lengthSpikes])-1] = spikes[j];
                    counter++;
                }

                if (counter == numCells)
                {
                    tempOutput[row] = spikes[j];
                    tempOutput[row+lengthSpikes] = spikes[i];
                    oldEnd = i;
                    row++;
                    break;
                }
                j--;
            }

            if (counter != numCells)            
            {
                /* Now look forwards for the missing pieces */
                j = i+1;
                while (((spikes[j] - spikes[i]) < windowSize) && (j < lengthSpikes) && (spikes[j+lengthSpikes] != minCellID))
                {
                    if (cellOccurrances[((int) spikes[j+lengthSpikes])-1] == 0)
                    {                        
                        counter++;
                    }
                    cellOccurrances[((int) spikes[j+lengthSpikes])-1] = spikes[j];

                    if (counter == numCells)
                    {
                        /* Ensure all values are within reach */
                        currMin = spikes[j];
                        for (k = 0; k < numCells; k++)
                        {
                            if (cellOccurrances[k] < currMin){currMin = cellOccurrances[k];}
                            if (spikes[j] - cellOccurrances[k] > windowSize)
                            {
                                counter--;
                                cellOccurrances[k] = 0;
                            }
                        }
                        if (counter == numCells)
                        {
                            tempOutput[row] = currMin;
                            tempOutput[row+lengthSpikes] = spikes[j];
                            oldEnd = j;
                            row++;
                            exitLoop = 1;
                            break;
                        }
                    }
                    j++;
                }     
            }
        }
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
