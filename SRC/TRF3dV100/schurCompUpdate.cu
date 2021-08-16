#include "superlu_ddefs.h"
#include "lupanels_GPU.cuh"
#include "lupanels.hpp"
// this should be a device code
__device__
int_t lpanelGPU_t::find(int_t k)
{
    return 0;
}


__device__
int_t upanelGPU_t::find(int_t k)
{
    return 0;
}


__device__
int computeIndirectMapGPU(int* rcS2D,  int_t srcLen, int_t *srcVec,
                                         int_t dstLen, int_t *dstVec,
                                         int *dstIdx)
{
    int threadId = threadIdx.x;
    if (dstVec == NULL) /*uncompressed dimension*/
    {
        if(threadId < srcLen)
            rcS2D[threadId] = srcVec[threadId];
        __syncthreads();
        return 0; 
    }
    
    
    if(threadId < dstLen)
        dstIdx[dstVec[threadId]] = threadId;
    __syncthreads();
    
    if(threadId < srcLen)
        rcS2D[threadId] = dstIdx[srcVec[threadId]];
    __syncthreads();
    
    return 0;
}

__global__
void scatterGPU(
    int iSt,  int jSt, 
    double* gemmBuff, int LDgemmBuff,
    lpanelGPU_t lpanel, upanelGPU_t upanel, 
    LUstructGPU_t* dA)
{

    // calculate gi,gj
    int ii = iSt + blockIdx.x; 
    int jj = jSt + blockIdx.y; 

    int gi = lpanel.gid(ii);
    int gj = upanel.gid(jj);

    double *Dst;
    int_t lddst;
    int_t dstRowLen, dstColLen;
    int_t *dstRowList;
    int_t *dstColList;

    if (gj > gi) // its in upanel
    {
        int li = dA->g2lRow(gi);
        int lj = dA->uPanelVec[li].find(gj);
        Dst = dA->uPanelVec[li].blkPtr(lj);
        lddst = dA->supersize(gi);
        dstRowLen = dA->supersize(gi);
        dstRowList = NULL;
        dstColLen = dA->uPanelVec[li].nbcol(lj);
        dstColList = dA->uPanelVec[li].colList(lj);
        // std::cout<<li<<" "<<lj<<" Dst[0] is"<<Dst[0] << "\n";
    }
    else
    {
        int lj = dA->g2lCol(gj);
        int li = dA->lPanelVec[lj].find(gi);
        Dst = dA->lPanelVec[lj].blkPtr(li);
        lddst = dA->lPanelVec[lj].LDA();
        dstRowLen = dA->lPanelVec[lj].nbrow(li);
        dstRowList = dA->lPanelVec[lj].rowList(li);
        dstColLen = dA->supersize(gj);
        dstColList = NULL;
    }

    // compute source row to dest row mapping
    int maxSuperSize = dA->maxSuperSize; 
    extern __shared__ int baseSharedPtr[]; 
    int* rowS2D = baseSharedPtr;
    int* colS2D = &rowS2D[maxSuperSize];
    int* dstIdx = &colS2D[maxSuperSize];

    int nrows = lpanel.nbrow(ii);
    int ncols = upanel.nbcol(jj);
    // lpanel.rowList(ii), upanel.colList(jj)
    
    computeIndirectMapGPU(rowS2D,  nrows, lpanel.rowList(ii),
        dstRowLen, dstRowList, dstIdx);
    
// compute source col to dest col mapping
    computeIndirectMapGPU(colS2D, ncols, upanel.colList(jj),
        dstColLen, dstColList, dstIdx);
    

    int nThreads = blockDim.x; 
    int colsPerThreadBlock = nThreads/ nrows;
    
    int threadId = threadIdx.x;
    


    double* Src = &gemmBuff[ lpanel.stRow(ii)+ 
            upanel.stCol(jj)* LDgemmBuff];
    int ldsrc = LDgemmBuff; 
    // TODO: this seems inefficient 
    if (threadId < nrows * colsPerThreadBlock)
	{
		/* 1D threads are logically arranged in 2D shape. */
		int i = threadId % nrows;
		int j = threadId / nrows;

        #pragma unroll 4
        while(j<ncols)
        {
            Dst[rowS2D[i] + lddst * colS2D[j]] -= Src[i + ldsrc * j];
            j += colsPerThreadBlock;
        }
		
	}

} 


int_t LUstruct_v100::dSchurComplementUpdate(
    cublasHandle_t handle, int streamId, 
    int_t k, lpanel_t &lpanel, upanel_t &upanel)
{
    // TODO: redefine isEmpty so this works out 
    if (lpanel.isEmpty() || upanel.isEmpty())
        return 0;

    int_t st_lb = 0;
    if (myrow == krow(k))
        st_lb = 1;

    int_t nlb = lpanel.nblocks();
    int_t nub = upanel.nblocks();

    int iSt =st_lb;
    int iEnd =0; 
    
    

    int nrows = lpanel.stRow(nlb) - lpanel.stRow(st_lb);
    int ncols = upanel.nzcols(); 

    int maxGemmRows = nrows;  
    int maxGemmCols = ncols;
    // entire gemm doesn't fit in gemm buffer
    if(nrows* ncols > A_gpu.gemmBufferSize)
    {
        int maxGemmOpSize = (int) sqrt(A_gpu.gemmBufferSize);
        int numberofRowChunks = (nrows +maxGemmOpSize-1)/maxGemmOpSize;
        maxGemmRows =   nrows /numberofRowChunks;
        maxGemmCols = A_gpu.gemmBufferSize/ maxGemmRows; 
    }
    
    while(iEnd< nlb)
    {
        iSt = iEnd;
        iEnd = lpanel.getEndBlock(iSt, maxGemmRows);
        
        
        int jSt =0;
        int jEnd =0; 
        while(jEnd< nub)
        {
            jSt = jEnd; 
            jEnd = upanel.getEndBlock(jSt, maxGemmCols);
            cudaStream_t cuStream = A_gpu.cuStreams[streamId];
            cublasSetStream(handle, cuStream);
            int gemm_m = lpanel.stRow(iEnd) - lpanel.stRow(iSt);
            int gemm_n = upanel.stCol(jEnd) - upanel.stCol(jSt);
            int gemm_k = supersize(k);
            double alpha = 1.0;
            double beta = 0.0; 
            
            cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                        gemm_m, gemm_n, gemm_k, &alpha,
                        lpanel.blkPtrGPU(iSt), lpanel.LDA(),
                        upanel.blkPtrGPU(jSt), upanel.LDA(), &beta,
                        A_gpu.gpuGemmBuffs[streamId], gemm_m);


            // setting up scatter 
            dim3 dimBlock(ldt); // 1d thread
            dim3 dimGrid(iEnd - iSt, jEnd - jSt);
            size_t sharedMemorySize=3* A_gpu.maxSuperSize; 

            scatterGPU<<<dimGrid, dimBlock, sharedMemorySize, cuStream>>>(
                iSt, jSt, 
                A_gpu.gpuGemmBuffs[streamId], gemm_m,
                lpanel.gpuPanel, upanel.gpuPanel, dA_gpu);   
            
            
		}
    }
    return 0;
}




