#include "superlu_ddefs.h"

sForest_t **compute_sForests(int_t nsupers,  dLUstruct_t *LUstruct, gridinfo3d_t *grid3d)
{
    // Calculation of supernodal etree
    int_t *setree = supernodal_etree(nsupers, LUstruct->etree, LUstruct->Glu_persist->supno, LUstruct->Glu_persist->xsup);

    // Conversion of supernodal etree to list
    treeList_t *treeList = setree2list(nsupers, setree);

    // Calculation of tree weight
    calcTreeWeight(nsupers, setree, treeList, LUstruct->Glu_persist->xsup);

    // Calculation of maximum level
    int_t maxLvl = log2i(grid3d->zscp.Np) + 1;

    // Generation of forests
    sForest_t **sForests = getForests(maxLvl, nsupers, setree, treeList);

    // Allocate trf3d data structure
    // LUstruct->trf3Dpart = (dtrf3Dpartition_t *)SUPERLU_MALLOC(sizeof(dtrf3Dpartition_t));
    // dtrf3Dpartition_t *trf3Dpart = LUstruct->trf3Dpart;

    return sForests;
}


gEtreeInfo_t fillEtreeInfo( int_t nsupers, int_t* setree, treeList_t *treeList) {
    
    gEtreeInfo_t gEtreeInfo;
    gEtreeInfo.setree = setree;
    gEtreeInfo.numChildLeft = (int_t *)SUPERLU_MALLOC(sizeof(int_t) * nsupers);
    for (int_t i = 0; i < nsupers; ++i)
    {
        gEtreeInfo.numChildLeft[i] = treeList[i].numChild;
    }

    return gEtreeInfo;
}

int_t* create_iperm_c_supno(int_t nsupers, superlu_dist_options_t *options, 
    dLUstruct_t *LUstruct, gridinfo3d_t *grid3d) {
    gridinfo_t *grid = &(grid3d->grid2d);
    int_t *perm_c_supno = getPerm_c_supno(nsupers, options,
                                          LUstruct->etree,
                                          LUstruct->Glu_persist,
                                          LUstruct->Llu->Lrowind_bc_ptr,
                                          LUstruct->Llu->Ufstnz_br_ptr, grid);
    int_t *iperm_c_supno = getFactIperm(perm_c_supno, nsupers);
    SUPERLU_FREE(perm_c_supno);
    return iperm_c_supno;
}

int_t *createSupernode2TreeMap(int_t nsupers, int_t maxLvl, int_t *gNodeCount, int_t **gNodeLists)
{
    int_t *supernode2treeMap = SUPERLU_MALLOC(nsupers * sizeof(int_t));
    int_t numForests = (1 << maxLvl) - 1;

    for (int_t Fr = 0; Fr < numForests; ++Fr)
    {
        for (int_t nd = 0; nd < gNodeCount[Fr]; ++nd)
        {
            supernode2treeMap[gNodeLists[Fr][nd]] = Fr;
        }
    }

    return supernode2treeMap;
}

SupernodeToGridMap_t* createSuperGridMap(int_t nsuper,int_t maxLvl, int_t *myTreeIdxs, 
    int_t *myZeroTrIdxs, int_t* gNodeCount, int_t** gNodeLists)
{
    SupernodeToGridMap_t* superGridMap = SUPERLU_MALLOC(nsuper * sizeof(SupernodeToGridMap_t));
    for (int_t i = 0; i < nsuper; ++i)
    {
        /* initialize with NOT_IN_GRID */
        superGridMap[i] = NOT_IN_GRID;
    }

    
    for(int_t lvl =0; lvl<maxLvl; lvl++ )
    {
        int_t treeIdx = myTreeIdxs[lvl];
        int_t zeroTrIdx = myZeroTrIdxs[lvl];

        for (int_t nd = 0; nd < gNodeCount[treeIdx]; ++nd)
        {
            /* code */
            if(zeroTrIdx)
                superGridMap[gNodeLists[treeIdx][nd]] = IN_GRID_ZERO;
            else
                superGridMap[gNodeLists[treeIdx][nd]] = IN_GRID_AIJ;
        }
    }
    return superGridMap;
}
void newTrfPartitionInit(int_t nsupers,  dLUstruct_t *LUstruct, gridinfo3d_t *grid3d)
{

    gridinfo_t* grid = &(grid3d->grid2d);
    int iam = grid3d->iam;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC (iam, "Enter newTrfPartitionInit()");
#endif

    // check parameters
    if (LUstruct->trf3Dpart == NULL || grid3d == NULL)
    {
        fprintf(stderr, "Error: Invalid arguments to newTrfPartitionInit().\n");
        return;
    }

      // Calculation of supernodal etree
    int_t *setree = supernodal_etree(nsupers, LUstruct->etree, LUstruct->Glu_persist->supno, LUstruct->Glu_persist->xsup);

    // Conversion of supernodal etree to list
    treeList_t *treeList = setree2list(nsupers, setree);

// YL: The essential difference between this function and dinitTrf3Dpartition_allgrid to avoid calling pddistribute* twice is that Piyush has removed the treelist weight update function below (and iperm_c_supno as well), which requires the LU data structure  
#if 0 
    /*update treelist with weight and depth*/
    getSCUweight_allgrid(nsupers, treeList, xsup,
        LUstruct->Llu->Lrowind_bc_ptr, LUstruct->Llu->Ufstnz_br_ptr,
        grid3d);
#endif
    // Calculation of tree weight
    calcTreeWeight(nsupers, setree, treeList, LUstruct->Glu_persist->xsup);

    // Calculation of maximum level
    int_t maxLvl = log2i(grid3d->zscp.Np) + 1;

    // Generation of forests
    sForest_t **sForests = getForests(maxLvl, nsupers, setree, treeList);

    dtrf3Dpartition_t *trf3Dpart = LUstruct->trf3Dpart;
    trf3Dpart->sForests = sForests;
    trf3Dpart->nsupers = nsupers;
      int_t *myTreeIdxs = getGridTrees(grid3d);
    int_t *myZeroTrIdxs = getReplicatedTrees(grid3d);
    int_t *gNodeCount = getNodeCountsFr(maxLvl, sForests);
    int_t **gNodeLists = getNodeListFr(maxLvl, sForests); // reuse NodeLists stored in sForests[]

    // dinit3DLUstructForest(myTreeIdxs, myZeroTrIdxs,
    //                       sForests, LUstruct, grid3d);
    int_t *myNodeCount = getMyNodeCountsFr(maxLvl, myTreeIdxs, sForests);
    int_t **treePerm = getTreePermFr(myTreeIdxs, sForests, grid3d);
    int* supernodeMask = SUPERLU_MALLOC(nsupers*sizeof(int));
    for (int ii = 0; ii < nsupers; ++ii)
        supernodeMask[ii]=0;
    for (int lvl = 0; lvl < maxLvl; ++lvl)
    {
        // printf("iam %5d lvl %5d myNodeCount[lvl] %5d\n",grid3d->iam, lvl,myNodeCount[lvl]);
        for (int nd = 0; nd < myNodeCount[lvl]; ++nd)
        {
            supernodeMask[treePerm[lvl][nd]]=1;
        }
    }





    // dLUValSubBuf_t *LUvsb = SUPERLU_MALLOC(sizeof(dLUValSubBuf_t));
    // dLluBufInit(LUvsb, LUstruct);

#if (DEBUGlevel>=1)    
    // let count sum of gnodecount
    int_t gNodeCountSum = 0;
    for (int_t i = 0; i < (1 << maxLvl) - 1; ++i)
    {
        gNodeCountSum += gNodeCount[i];
    }
    printf(" Iam: %d, Nsupers %d, gnodecountSum =%d \n", grid3d->iam, nsupers, gNodeCountSum);
#endif

    /* Sherry 2/17/23
       Compute buffer sizes needed for diagonal LU blocks and C matrices in GEMM. */

    
    iam = grid->iam;  /* 'grid' is 2D grid */
    int k, k0, k_st, k_end, offset, nsupc, krow, kcol;
    int myrow = MYROW (iam, grid);
    int mycol = MYCOL (iam, grid);
    int_t *xsup  = LUstruct->Glu_persist->xsup;
    
#if 0    
    int krow = PROW (k, grid);
    int kcol = PCOL (k, grid);
    int_t** Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
    double** Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;

    int_t** Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
    double** Unzval_br_ptr = Llu->Unzval_br_ptr;
#endif    
    
    int mxLeafNode = 0; // Yang: only need to check the leaf level of topoInfo as the factorization proceeds level by level 
    for (int ilvl = 0; ilvl < maxLvl; ++ilvl) {
        if (sForests[myTreeIdxs[ilvl]] && sForests[myTreeIdxs[ilvl]]->topoInfo.eTreeTopLims[1] > mxLeafNode )
            mxLeafNode    = sForests[myTreeIdxs[ilvl]]->topoInfo.eTreeTopLims[1];
    }

    // Yang: use ldts to track the maximum needed buffer sizes per node of topoInfo 
    //int *ldts = (int*) SUPERLU_MALLOC(mxLeafNode*sizeof(int));
    //for (int i = 0; i < mxLeafNode; ++i) {  //????????
    //ldts[i]=1;
    //}
    int *ldts = int32Calloc_dist(mxLeafNode);

    for (int ilvl = 0; ilvl < maxLvl; ++ilvl) {  /* Loop through the Pz tree levels */
        int treeId = myTreeIdxs[ilvl];
        sForest_t* sforest = sForests[treeId];
        if (sforest){
            int_t *perm_node = sforest->nodeList ; /* permuted list, in order of factorization */
	    int maxTopoLevel = sforest->topoInfo.numLvl;/* number of levels at each outer-tree node */
            for (int topoLvl = 0; topoLvl < maxTopoLevel; ++topoLvl)
            {
                /* code */
                k_st = sforest->topoInfo.eTreeTopLims[topoLvl];
                k_end = sforest->topoInfo.eTreeTopLims[topoLvl + 1];
		//printf("\t..topoLvl %d, k_st %d, k_end %d\n", topoLvl, k_st, k_end);
		
                for (int k0 = k_st; k0 < k_end; ++k0)
                {
                    offset = k0 - k_st;
                    k = perm_node[k0];
                    nsupc = (xsup[k+1]-xsup[k]);
                    krow = PROW (k, grid);
                    kcol = PCOL (k, grid);
                    if ( myrow == krow || mycol == kcol )  /* diagonal process */
                    {
		        ldts[offset] = SUPERLU_MAX(ldts[offset], nsupc);
                    }
#if 0 /* GPU gemm buffers can only be set on GPU side, because here we only know
	 the size of U data structure on CPU.  It is different on GPU */
                    if ( mycol == kcol ) { /* processes owning L panel */
		      
		    }
                    if ( myrow == krow ) 
			gemmCsizes[offset] = SUPERLU_MAX(ldts[offset], ???);
#endif		    
                }               
            }
        }
    }




    trf3Dpart->gEtreeInfo = fillEtreeInfo(nsupers, setree, treeList);
    // trf3Dpart->iperm_c_supno = iperm_c_supno;
    trf3Dpart->myNodeCount = myNodeCount;
    trf3Dpart->myTreeIdxs = myTreeIdxs;
    trf3Dpart->myZeroTrIdxs = myZeroTrIdxs;
    trf3Dpart->sForests = sForests;
    trf3Dpart->treePerm = treePerm;
    trf3Dpart->maxLvl = maxLvl;
    // trf3Dpart->LUvsb = LUvsb;
    trf3Dpart->supernode2treeMap = createSupernode2TreeMap(nsupers, maxLvl, gNodeCount, gNodeLists);
    trf3Dpart->superGridMap = createSuperGridMap(nsupers, maxLvl, myTreeIdxs, myZeroTrIdxs, gNodeCount, gNodeLists);
    trf3Dpart->supernodeMask = supernodeMask;
    trf3Dpart->mxLeafNode = mxLeafNode;  // Sherry added these 3
    trf3Dpart->diagDims = ldts;
    //trf3Dpart->gemmCsizes = gemmCsizes;
    // Sherry added
    // Deallocate storage
    SUPERLU_FREE(gNodeCount);
    SUPERLU_FREE(gNodeLists);
    free_treelist(nsupers, treeList);
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC (iam, "Exit newTrfPartitionInit()");
#endif

}

/**
 * @file allocBcastArray.c
 * @brief Function to allocate and broadcast an array in a MPI environment.
 */

#include <mpi.h>
#include <stdlib.h>

/**
 * @brief Allocates and broadcasts an array in a MPI environment.
 *
 * This function sends the size from the root process to all other processes in
 * the communicator. If the process is not the root, it receives the size from
 * the root and allocates the array. Then, the function broadcasts the array
 * from the root process to all other processes in the communicator.
 *
 * @param array Pointer to the array to be allocated and broadcasted.
 * @param size The size of the array.
 * @param comm The MPI communicator.
 * @param root The root process.
 */
void allocBcastArray(void **array, int size, int root, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank); // Get the rank of the current process

    // Check if the size is valid
    if (rank == root)
    {
        if (size <= 0)
        {
            fprintf(stderr, "Error: Size should be a positive integer.\n");
            MPI_Abort(comm, EXIT_FAILURE);
        }
        
    }

    // Send the size from root to all other processes in the communicator
    MPI_Bcast(&size, 1, MPI_INT, root, comm);

    // If I am not the root, receive the size from the root and allocate the array
    if (rank != root)
    {
        *array = SUPERLU_MALLOC(size);
        if (*array == NULL)
        {
            fprintf(stderr, "Error: Failed to allocate memory.\n");
            MPI_Abort(comm, EXIT_FAILURE);
        }
    }

    // Then broadcast the array from the root to all other processes
    MPI_Bcast(*array, size, MPI_BYTE, root, comm);
}

// function to broad permuted sparse matrix and symbolic factorization data from
// 2d to 3d grid

void bcastPermutedSparseA(SuperMatrix *A, 
                          dScalePermstruct_t *ScalePermstruct,
                          Glu_freeable_t *Glu_freeable, 
                          dLUstruct_t *LUstruct, gridinfo3d_t *grid3d)
{
    int_t m = A->nrow;
	int_t n = A->ncol;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    dLocalLU_t *Llu = LUstruct->Llu;
    NRformat_loc *Astore   = (NRformat_loc *) A->Store;
    // check if the varaibles are not NULL
    if (A == NULL || ScalePermstruct == NULL ||
        Glu_freeable == NULL || LUstruct == NULL  || grid3d == NULL ||  
        Glu_persist == NULL || Llu == NULL || Astore == NULL)
    {
        fprintf(stderr, "Error: Invalid arguments to bcastPermutedSparseA().\n");
        return;
    }

    /* broadcast etree */
    int_t *etree = LUstruct->etree;
    MPI_Bcast( etree, n, mpi_int_t, 0,  grid3d->zscp.comm);
    

    // list of all the arrays to be broadcasted 
    // A, ScalePermstruct, Glu_freeable, LUstruct
    int_t nsupers;
    
    if (!grid3d->zscp.Iam)
        nsupers = Glu_persist->supno[n-1] + 1;
    // broadcast the number of supernodes
    MPI_Bcast(&nsupers, 1, mpi_int_t, 0, grid3d->zscp.comm);
    
    /* ==== Broadcasting GLU_persist   ======= */
    // what is the size of xsup and supno?
    allocBcastArray( &(Glu_persist->xsup), (nsupers+1)*sizeof(int_t), 
        0, grid3d->zscp.comm);
    allocBcastArray( &(Glu_persist->supno), (n)*sizeof(int_t), 
        0, grid3d->zscp.comm);
    int_t *xsup = Glu_persist->xsup;    /* supernode and column mapping */
    int_t *supno = Glu_persist->supno;
    

    /* ==== Broadcasting ScalePermstruct ======= */
//     typedef struct {
//     DiagScale_t DiagScale; // enum 1
//     double *R;  (double*) dimension (A->nrow)
//     double *C;   (double*) dimension (A->ncol)
//     int_t  *perm_r; (int_t*) dimension (A->nrow)
//     int_t  *perm_c; (int_t*) dimension (A->ncol)
// } dScalePermstruct_t;

    MPI_Bcast(&(ScalePermstruct->DiagScale), sizeof(DiagScale_t), MPI_BYTE, 0, grid3d->zscp.comm);
    
/***** YL: remove the allocation in the following as perm_r/perm_c has been allocated on all grids by dScalePermstructInit
*/
#if 1  
    MPI_Bcast(ScalePermstruct->perm_r, m*sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);  
    MPI_Bcast(ScalePermstruct->perm_c, n*sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);  
#else
    allocBcastArray ( &(ScalePermstruct->perm_r), m*sizeof(int_t), 
        0, grid3d->zscp.comm);
    allocBcastArray ( &(ScalePermstruct->perm_c), n*sizeof(int_t),
        0, grid3d->zscp.comm);
#endif

    allocBcastArray ( &(ScalePermstruct->R), m*sizeof(double),
        0, grid3d->zscp.comm);
    allocBcastArray ( &(ScalePermstruct->C), n*sizeof(double),
        0, grid3d->zscp.comm);
    

    /* ==== Broadcasting Glu_freeable ======= */
//     typedef struct {
//     int_t     *lsub;     /* compressed L subscripts */
//     int_t     *xlsub;        // i think its size is n+1
//     int_t     *usub;     /* compressed U subscripts i think nzumax*/
//     int_t     *xusub;
//     int_t     nzlmax;    /* current max size of lsub */
//     int_t     nzumax;    /*    "    "    "      usub */
//     LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
//     //int_t     *llvl;     /* keep track of level in L for level-based ILU */
//     //int_t     *ulvl;     /* keep track of level in U for level-based ILU */
//     int64_t nnzLU;   /* number of nonzeros in L+U*/
// } Glu_freeable_t;

    allocBcastArray( &(Glu_freeable->lsub), Glu_freeable->nzlmax*sizeof(int_t),
        0, grid3d->zscp.comm);
    allocBcastArray( &(Glu_freeable->xlsub), (n+1)*sizeof(int_t),
        0, grid3d->zscp.comm);
    allocBcastArray( &(Glu_freeable->usub), Glu_freeable->nzumax*sizeof(int_t),
        0, grid3d->zscp.comm);
    allocBcastArray( &(Glu_freeable->xusub), (n+1)*sizeof(int_t),
        0, grid3d->zscp.comm);
    MPI_Bcast(&(Glu_freeable->nzlmax), sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(&(Glu_freeable->nzumax), sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(&(Glu_freeable->nnzLU), sizeof(int64_t), MPI_BYTE, 0, grid3d->zscp.comm);

    if(grid3d->zscp.Iam)
    {
        Glu_freeable->MemModel = SYSTEM;
    }
    
    
    /* ==== Broadcasting permuted sparse matrix ======= */
    // Astore = (NRformat_loc *) A->Store;
//     /typedef struct {
//     int_t nnz_loc;   /* number of nonzeros in the local submatrix */
//     int_t m_loc;     /* number of rows local to this processor */
//     int_t fst_row;   /* global index of the first row */
//     void  *nzval;    /* pointer to array of nonzero values, packed by row */
//     int_t *rowptr;   /* pointer to array of beginning of rows in nzval[] 
// 			and colind[]  */
//     int_t *colind;   /* pointer to array of column indices of the nonzeros */
//                      /* Note:
// 			Zero-based indexing is used;
// 			rowptr[] has n_loc + 1 entries, the last one pointing
// 			beyond the last row, so that rowptr[n_loc] = nnz_loc.*/
// } NRformat_loc;

    
    // NRformat_loc *Astore = (NRformat_loc *) A->Store;
    MPI_Bcast(&(Astore->nnz_loc), sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(&(Astore->m_loc), sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(&(Astore->fst_row), sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);


/***** YL: remove the allocation in the following as dGatherNRformat_loc3d_allgrid instead of dGatherNRformat_loc3d has been called, which already allocate A->Store on all grids
 * Note the the broadcast is still needed as the A->Store has been scaled by scaleMatrixDiagonally only on grid 0
*/
#if 1
    MPI_Bcast(Astore->nzval, Astore->nnz_loc*sizeof(double), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(Astore->rowptr, (Astore->m_loc+1)*sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
    MPI_Bcast(Astore->colind, Astore->nnz_loc*sizeof(int_t), MPI_BYTE, 0, grid3d->zscp.comm);
#else    
    allocBcastArray( &(Astore->nzval), Astore->nnz_loc*sizeof(double),
        0, grid3d->zscp.comm);
    allocBcastArray( &(Astore->rowptr), (Astore->m_loc+1)*sizeof(int_t), 
        0, grid3d->zscp.comm);
    allocBcastArray( &(Astore->colind), Astore->nnz_loc*sizeof(int_t), 
        0, grid3d->zscp.comm);
#endif

}


