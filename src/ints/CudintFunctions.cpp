#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CudintF.h"

using namespace std;

struct Node {
  int index_a;
  int index_b;
  int index_r;
  int index_s;
  Node* next;
};

// only for the 1st Node
void initNode(struct Node *head,int a, int b, int r, int s){
  head->index_a = a;
  head->index_b = b;
  head->index_r = r;
  head->index_s = s;
  head->next =NULL;
}

// apending
void addNode(struct Node *head, int a, int b, int r, int s) {
  Node *newNode = new Node;
  newNode->index_a = a;
  newNode->index_b = b;
  newNode->index_r = r;
  newNode->index_s = s;
  newNode->next = NULL;

  Node *cur = head;
  while(cur) {
    if(cur->next == NULL) {
      cur->next = newNode;
      return;
    }
    cur = cur->next;
  }
}

void display(struct Node *head) {
  Node *list = head;
  printf("Contraction indices:\n");
  while(list) {
    printf("(%d,%d|%d,%d)\n", list->index_a, 
	   list->index_b, list->index_r, list->index_s);
    list = list->next;
  }
  printf("****************\n");
}


struct primitives {
  int contID;
  int index_i;
  int index_j;
  int index_k;
  int index_l;
  int cart_i;
  int cart_j;
  int cart_k;
  int cart_l;
  primitives* next;
};

// only for the 1st Node
void initPrimitives(struct primitives *primHead, int contID, int i, int j, int k, int l, int cart_i, int cart_j, int cart_k, int cart_l){
  primHead->contID = contID;
  primHead->index_i = i;
  primHead->index_j = j;
  primHead->index_k = k;
  primHead->index_l = l;
  primHead->cart_i = cart_i;
  primHead->cart_j = cart_j;
  primHead->cart_k = cart_k;
  primHead->cart_l = cart_l;
  primHead->next =NULL;
}

// apending
void addPrimitives(struct primitives *primHead, int contID, int i, int j, int k, int l, int cart_i, int cart_j, int cart_k, int cart_l) {
  primitives *newPrimitives = new primitives;
  newPrimitives->contID = contID;
  newPrimitives->index_i = i;
  newPrimitives->index_j = j;
  newPrimitives->index_k = k;
  newPrimitives->index_l = l;
  newPrimitives->cart_i = cart_i;
  newPrimitives->cart_j = cart_j;
  newPrimitives->cart_k = cart_k;
  newPrimitives->cart_l = cart_l;
  newPrimitives->next =NULL;

  primitives *cur = primHead;
  while(cur) {
    if(cur->next == NULL) {
      cur->next = newPrimitives;
      return;
    }
    cur = cur->next;
  }
}

void displayPrimitives(struct primitives *primHead) {
  primitives *list = primHead;
  printf("Primitives indices:\n");
  while(list) {
    printf("[%d]: (%d,%d|%d,%d) -> cartesians: [%d,%d|%d,%d]\n", 
	   list->contID,
	   list->index_i, 
	   list->index_j, 
	   list->index_k,
	   list->index_l,
	   list->cart_i,
	   list->cart_j,
	   list->cart_k,
	   list->cart_l);
    list = list->next;
  }
  printf("****************\n");
}

void getIndices(int numberOfContractions, int *contLength, int *angularMoments, int *numCartesianOrbitals, int *labelsForContractions, int *totalPrimitives)
{
  int m, auxCount, auxCount2, auxCount3;
  int a, b, r, s, n, u;
  int aa, bb, rr, ss;
  int i,j,k,l, counterPerCatersian;
  int ii, jj, kk, ll;
  int pa, pb, pr, ps;
  int * pi, * pj, * pk, * pl;
  int * poi, * poj, * pok, * pol;
  int unicintegrals, unicintegralsMem;
  int *numberOfPPUC, *contIndices, *realContIndices, *primIndices;
  int contIndicesMem, primIndicesMem;
  int indexA, indexB, indexR, indexS;
  int *auxPrimitiveptr = totalPrimitives;
  int primitiveCounter;

  struct Node *newHead;
  struct Node *head = new Node;
  struct Node *newHeadB;
  struct Node *headB = new Node;
  struct Node *newHeadC;
  struct Node *headC = new Node;

  struct primitives *newprimHead;
  struct primitives *primHead = new primitives;

  // unicintegrals = ((numberOfContractions*(numberOfContractions+1)/2)+1)*(numberOfContractions*(numberOfContractions+1)/2)/2;

  // unicintegralsMem = unicintegrals*sizeof(int);
  auxCount2=0;
  auxCount3=0;
  m=0;
  counterPerCatersian = 0;
  *auxPrimitiveptr = 0;
  primitiveCounter = 0;
  for( a = 1;  a<=numberOfContractions; a++)
    {
      n = a;
      for( b = a; b<=numberOfContractions;b++)
  	{
          u = b;
          for( r = n ;r <=numberOfContractions;r++)
  	    {
  	      for( s = u; s<=numberOfContractions; s++)
  		{

		  *auxPrimitiveptr += contLength[a-1]*contLength[b-1]*contLength[r-1]*contLength[s-1];
						 
  		  int aux = 0;
  		  int order = 0;

  		  //permuted index
  		  aa = a;
  		  bb = b;
  		  rr = r;
  		  ss = s;
                   
                  //pointer to permuted index under a not permuted loop
  		  pi = &i;
  		  pj = &j;
  		  pk = &k;
  		  pl = &l;
                   
  		  //pointer to not permuted index under a permuted loop
  		  poi = &i;
  		  poj = &j;
  		  pok = &k;
  		  pol = &l;
		  
  		  if (angularMoments[a-1] < angularMoments[b-1])
  		    {
  		      aa = b;
  		      bb = a;
                      
  		      pi = &j;
  		      pj = &i;
                      
  		      poi = &j;
  		      poj = &i;
                      
  		      order = order + 1;
  		    }
                   
  		  if (angularMoments[r-1] < angularMoments[s-1])
  		    {
  		      // printf("encontre un r menor: %d , %d\n", angularMoments[r-1], angularMoments[s-1]);
  		      rr = s;
  		      ss = r;
                      
  		      pk = &l;
  		      pl = &k;
                      
  		      pok = &l;
  		      pol = &k;
                      
  		      order = order + 3;
  		    }
                   
                   if((angularMoments[a-1] + angularMoments[b-1]) < (angularMoments[r-1] + angularMoments[s-1]))
  		     {  
  		       aux = aa;
  		       aa = rr;
  		       rr = aux;
                      
  		       aux = bb;
  		       bb = ss;
  		       ss = aux;
                      
  		       switch(order)
  			 {
  			 case 0:
  			   pi = &k;
  			   pj = &l;
  			   pk = &i;
  			   pl = &j;
                         
  			   poi = &k;
  			   poj = &l;
  			   pok = &i;
  			   pol = &j;
  			   break;
  			 case 1:
  			   pi = &k;
  			   pj = &l;
  			   pk = &j;
  			   pl = &i;
                         
  			   poi = &l;
  			   poj = &k;
  			   pok = &i;
  			   pol = &j;
  			   break;
  			 case 3:
  			   pi = &l;
  			   pj = &k;
  			   pk = &i;
  			   pl = &j;
                         
  			   poi = &k;
  			   poj = &l;
  			   pok = &j;
  			   pol = &i;
  			   break;
  			 case 4:
  			   pi = &l;
  			   pj = &k;
  			   pk = &j;
  			   pl = &i;
                         
  			   poi = &l;
  			   poj = &k;
  			   pok = &j;
  			   pol = &i;
  			   break;
  			 }
  		     }

		   for(i=1; i<=numCartesianOrbitals[aa-1]; i++)
		     for(j=1; j<=numCartesianOrbitals[bb-1]; j++)
		       for(k=1; k<=numCartesianOrbitals[rr-1]; k++)
			 for(l=1; l<=numCartesianOrbitals[ss-1]; l++)
			   {

			     // index not permuted
			     pa=labelsForContractions[a-1] + *poi - 1;
			     pb=labelsForContractions[b-1] + *poj - 1;
			     pr=labelsForContractions[r-1] + *pok - 1;
			     ps=labelsForContractions[s-1] + *pol - 1;

                               
			     if( pa <= pb && pr <= ps && (pa*1000)+pb >= (pr*1000)+ps)
			       {
				 aux = pa;
				 pa = pr;
				 pr = aux;
                                  
				 aux = pb;
				 pb = ps;
				 ps = aux;
			       }

			     if( pa <= pb && pr <= ps )
			       {
				 auxCount = 0;
				 auxCount = contLength[aa-1]*contLength[bb-1]*contLength[rr-1]*contLength[ss-1];
				 auxCount2 = auxCount3;
				 auxCount3 += auxCount;
				 if(counterPerCatersian==0)
				   {
				     initNode(head,aa,bb,rr,ss);
				     initNode(headB,pa,pb,pr,ps);
				     initNode(headC,auxCount, auxCount2, auxCount3, 0);
				   }
				 else
				   {
				     addNode(head,aa,bb,rr,ss);
				     addNode(headB,pa,pb,pr,ps);
				     addNode(headC,auxCount, auxCount2, auxCount3, 0);
				   }


				 for(ii=1;ii<=contLength[aa-1];ii++)
				   for(jj=1;jj<=contLength[bb-1];jj++)
				     for(kk=1;kk<=contLength[rr-1];kk++)
				       for(ll=1;ll<=contLength[ss-1];ll++)
					 {
					   if(primitiveCounter==0)
					     {
					       initPrimitives(primHead, 
							      counterPerCatersian,
							      ii, jj, kk, ll, 
							      i, j, k, l);
					     }
					   else
					     {
					       addPrimitives(primHead, 
							     counterPerCatersian,
							     ii, jj, kk, ll, 
							     i, j, k, l);
					     }
					   primitiveCounter++;
					 }

				 // printf("Contraction (%d) -> [%d]: (%d,%d|%d,%d) -> [%d, %d, %d, %d] (%d,%d|%d,%d) \n", 
				 // 	m, counterPerCatersian, 
				 // 	aa, bb, rr, ss, 
				 // 	i, j, k, l, 
				 // 	pa, pb, pr, ps);
				 counterPerCatersian++;
			       }
			   }
		   m++;	   
  		}
  	      u = r+1;
  	    }
  	}
    }

  // display(head);
  // display(headB);
  // displayPrimitives(primHead);

  primIndicesMem = primitiveCounter*sizeof(int);
  contIndicesMem = counterPerCatersian*sizeof(int);

  contIndices = (int *)malloc(4*contIndicesMem);
  realContIndices = (int *)malloc(4*contIndicesMem);
  //numberOfPPC = Number of Primitives per Unic Integral Contraction
  numberOfPPUC = (int *)malloc(3*contIndicesMem);

  primIndices = (int *)malloc(9*primIndicesMem);
  
  // printf("Total primitives new: %d\n", primitiveCounter);
  // printf("Contraction indices fake and real\n");
  i=0;
  while(head) {
    contIndices[i*4] = head->index_a;
    contIndices[i*4+1] = head->index_b;
    contIndices[i*4+2] = head->index_r;
    contIndices[i*4+3] = head->index_s;
    realContIndices[i*4] = headB->index_a;
    realContIndices[i*4+1] = headB->index_b;
    realContIndices[i*4+2] = headB->index_r;
    realContIndices[i*4+3] = headB->index_s;
    numberOfPPUC[i*3] = headC->index_a;
    numberOfPPUC[i*3+1] = headC->index_b;
    numberOfPPUC[i*3+2] = headC->index_r;
    // printf("[%d: %d->%d] (%d,%d|%d,%d) -> (%d,%d|%d,%d)\n", 
    // 	   numberOfPPUC[i*3],
    // 	   numberOfPPUC[i*3+1],
    // 	   numberOfPPUC[i*3+2],
    // 	   contIndices[i*4],
    // 	   contIndices[i*4+1],
    // 	   contIndices[i*4+2], 
    // 	   contIndices[i*4+3],
    // 	   realContIndices[i*4],
    // 	   realContIndices[i*4+1],
    // 	   realContIndices[i*4+2],
    // 	   realContIndices[i*4+3]);
    head = head->next;
    headB = headB->next;
    headC = headC->next;
    i++;
  }
  // printf("****************\n");  

  // printf("Primitive indices\n");
  i=0;
  while(primHead) {
    primIndices[i*9] = primHead->contID;
    primIndices[i*9+1] = primHead->index_i;
    primIndices[i*9+2] = primHead->index_j;
    primIndices[i*9+3] = primHead->index_k;
    primIndices[i*9+4] = primHead->index_l;
    primIndices[i*9+5] = primHead->cart_i ;
    primIndices[i*9+6] = primHead->cart_j ;
    primIndices[i*9+7] = primHead->cart_k ;
    primIndices[i*9+8] = primHead->cart_l ;
    // printf("[%d]: (%d,%d|%d,%d) -> [%d,%d|%d,%d]\n", 
    // 	   primIndices[i*9],
    // 	   primIndices[i*9+1],
    // 	   primIndices[i*9+2],
    // 	   primIndices[i*9+3],
    // 	   primIndices[i*9+4],
    // 	   primIndices[i*9+5],
    // 	   primIndices[i*9+6],
    // 	   primIndices[i*9+7],
    // 	   primIndices[i*9+8] );
    primHead = primHead->next;
    i++;
  }
  // printf("****************\n");  


}

