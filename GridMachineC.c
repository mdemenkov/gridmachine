
#include <math.h>
#include <stdio.h>
#include <mex.h>


typedef struct RealMatrixTag
{
	  int m,n,num,dis;
	  double *R;
}
  Matrix;

typedef struct BoolMatrixTag
{
	  int m,n,num,dis;
	  unsigned char *R;
}
  BoolMatrix;

typedef double (*UserFunc)(Matrix *x,int Component,double *FuncData);

struct InputStruc
{
  int SpaceDim,JacStep,Density,Progress,
	  SolSize,MaxTime;
  Matrix Limits,DeriSteps;
  double NewtonZero,SamePoints,Radius;
  UserFunc Func;
  double *FuncData;

}
  Input;

struct WorkTag
{
  int NVertices;
  long NCubes;

  BoolMatrix Vertices,Cube;
  Matrix XGrid,CubeVert,Range,CubeCenter,Sol,x;

}
  Work;

struct OutputStruc
{
  Matrix Sol,
	     SolutionsNumber,AllCubes,ProcessedCubes,NewtonFailed,RedundantSolutions;
}
  Output;

void ShowMat(Matrix *mat,char *str)
{ 
  int i,j,brow,pos;

  brow=0;mexPrintf("\n");mexPrintf(str);mexPrintf("\n");
  for(i=0;i<mat->m;i++)
  {
    pos=brow;
    for(j=0;j<mat->n;j++){mexPrintf("%4.3f  ",mat->R[pos]);pos+=mat->m;}
    brow+=1;mexPrintf("\n");
  }
  
}

void ShowNum(double val,char *str)
{
  mexPrintf(str);mexPrintf("%4.3f\n",val);
}

void ShowBoolMat(BoolMatrix *mat,char *str)
{ 
  int i,j,brow,pos;

  brow=0;mexPrintf("\n");mexPrintf(str);mexPrintf("\n");
  for(i=0;i<mat->m;i++)
  {
    pos=brow;
    for(j=0;j<mat->n;j++){mexPrintf("%1u ",mat->R[pos]);pos+=mat->m;}
    brow+=1;mexPrintf("\n");
  }
  
}

/*------------Function under investigation --------------------------*/

double example(Matrix *x,int Component,double *Data)
{
   double x1,x2,r;

   x1=x->R[x->dis+0];x2=x->R[x->dis+1];r=(*Data);

   
   switch(Component)
   {
      case 1:  return  x1*x1+x2*x2-r;
	  case 2:  return  x2;
   }
}

double bif3branch_onsphere(Matrix *x,int Component,double *Data)
{

	double x1,x2,x3,r;

	x1=x->R[x->dis+0];x2=x->R[x->dis+1];x3=x->R[x->dis+2];r=(*Data);

	switch(Component)
    {
   	  case 1: return x1*(x1*x1+x3);
      case 2: return x2*(x2*x2-x3);
      case 3: return x1*x1+x2*x2+x3*x3-r;
	}
}

double kubichek(Matrix *x,int Component,double *Data)
{
	double A,B,R0,D1,D2,cp;
	double X1,X2,X3,X4,Y1,Y2,Y3,Y4;

	A=2;B=6;R0=0.1;

	cp=(*Data);

	X1=x->R[x->dis+0];  Y1=x->R[x->dis+1];  X2=x->R[x->dis+2];  Y2=x->R[x->dis+3];
    X3=x->R[x->dis+4];  Y3=x->R[x->dis+5];  X4=x->R[x->dis+6];  Y4=x->R[x->dis+7];

    D1=pow(10,cp); 	D2=D1/R0;

    switch(Component)
	{
		case 1: return A-(B+1)*X1+X1*X1*Y1+D1*(X2-X1);
    
        case 2: return B*X1-X1*X1*Y1+D2*(Y2-Y1);
    
        case 3: return A-(B+1)*X2+X2*X2*Y2+D1*(X1-2*X2+X3);
    
        case 4: return B*X2-X2*X2*Y2+D2*(Y1-2*Y2+Y3);
    
        case 5: return A-(B+1)*X3+X3*X3*Y3+D1*(X2-2*X3+X4);
    
        case 6: return B*X3-X3*X3*Y3+D2*(Y2-2*Y3+Y4);
    
        case 7: return A-(B+1)*X4+X4*X4*Y4+D1*(X3-X4);
    
        case 8: return B*X4-X4*X4*Y4+D2*(Y3-Y4);
    }

}

/*------------ End function ------------------------------------------*/


/*-------------- Linear system solver from LAPACK -------------------------------*/

struct DGESV_DATA
{
  int N,NRHS,LDA,*IPIV,LDB,INFO;
  double *A,*B;
}
  dgesv_data;

extern void dgesv(int *N,int *NRHS,double *A,int *LDA,int *IPIV,double *B,
				  int *LDB,int *INFO);

void dgesv_init(int n)
{
	dgesv_data.N=n;
	dgesv_data.NRHS=1;
    dgesv_data.IPIV=(int *)mxCalloc(n,sizeof(int));
    dgesv_data.A=(double *)mxCalloc(n*n,sizeof(double));
	dgesv_data.LDA=n;
    dgesv_data.B=(double *)mxCalloc(n,sizeof(double));
    dgesv_data.LDB=n;
}

void dgesv_destroy()
{
	 mxFree(dgesv_data.IPIV);
	 mxFree(dgesv_data.A);
     mxFree(dgesv_data.B);
}

int SolveLinearSystem(Matrix *X,Matrix *A,Matrix *B)
{
     int i;

     for(i=0;i<A->num;i++) dgesv_data.A[i]=A->R[i];
     for(i=0;i<B->num;i++) dgesv_data.B[i]=B->R[i];
     dgesv(&dgesv_data.N,&dgesv_data.NRHS,dgesv_data.A,&dgesv_data.LDA,
		   dgesv_data.IPIV,dgesv_data.B,&dgesv_data.LDB,&dgesv_data.INFO);
	 if(dgesv_data.INFO==0)
	 {
        for(i=0;i<X->num;i++) X->R[i]=dgesv_data.B[i];
		return 1;
	 }
     else return 0;
}

/*----------------------End LAPACK ----------------------------------------------*/

void AllocMatrix(Matrix *mat,int m,int n)
{
	int num;

	mat->m=m;
	mat->n=n;
	num=m*n;
	mat->num=num;mat->dis=0;
	mat->R=(double *)mxCalloc(num,sizeof(double));
}

void AllocBoolMatrix(BoolMatrix *mat,int m,int n)
{
	int num;

	mat->m=m;
	mat->n=n;
	num=m*n;
	mat->num=num;mat->dis=0;
	mat->R=(unsigned char *)mxCalloc(num,sizeof(char));
}

void DestroyMatrix(Matrix *mat)
{
	mxFree(mat->R);
}

void DestroyBoolMatrix(BoolMatrix *mat)
{
	mxFree(mat->R);
}

void AddColumn(Matrix *out,Matrix *in)
{
  int i,j;

  for(j=out->dis,i=0;i<in->num;j++,i++) out->R[j]=in->R[i];
}

void CopyMat(Matrix *out,Matrix *in)
{
  int i;

  for(i=0;i<in->num;i++) out->R[i]=in->R[i];
}

void MatByNum(Matrix *mat,double scale)
{
  int i;

  for(i=0;i<mat->num;i++) mat->R[i]*=scale;
}

void MatByVec(Matrix *out,Matrix *mat,Matrix *in)
{
  int i,j,brow,pos;

  brow=0;
  for(i=0;i<mat->m;i++)
  {
    pos=brow;out->R[i]=0;
    for(j=0;j<mat->n;j++)
	{
		out->R[i]+=mat->R[pos]*in->R[j];
		pos+=mat->m;
	}
    brow+=1;
  }
}

void SubMat(Matrix *out,Matrix *in)
{
  int i;

  for(i=0;i<in->num;i++) out->R[i]-=in->R[i];
}

void Field2Matrix(mxArray *ptr,Matrix *mat,char *name)
{
	mxArray *field;

	field=mxGetField(ptr,0,name);
	mat->R=mxGetPr(field);
	mat->m=mxGetM(field);
	mat->n=mxGetN(field);
	mat->num=mat->n*mat->m;
}

double Field2Double(mxArray *ptr,char *name)
{
	mxArray *field;

	field=mxGetField(ptr,0,name);
	return (*mxGetPr(field));
}

void CreateOutField(mxArray *ptr,Matrix *mat,int m,int n,char *name)
{ 
     mxArray *field;

	 field=mxCreateDoubleMatrix(m,n,mxREAL);
	 mat->R=mxGetPr(field);
	 mat->m=m;
	 mat->n=n;
	 mat->num=n*m;
     mxSetField(ptr,0,name,field);
}


double InfNorm(Matrix *x)
{
  double mx,mdl;
  int i;

  mx=0;
  for(i=0;i<x->num;i++) {mdl=fabs(x->R[i]);if(mdl>mx) mx=mdl;}
  return mx;
}

int Contain(Matrix *mat,Matrix *x,double Tol)
{
    int flag,i,j;

	flag=1;
	for(i=mat->dis,j=0;j<x->m;i++,j++) 
		if(fabs(mat->R[i]-x->R[j])>Tol) {flag=0;break;}

	return flag;
}

double Min(Matrix *x)
{
  double mx;
  int i;

  mx=x->R[0];
  for(i=1;i<x->num;i++) if(mx>x->R[i]) mx=x->R[i];
  return mx;
}

double Max(Matrix *x)
{
  double mx;
  int i;

  mx=x->R[0];
  for(i=1;i<x->num;i++) if(mx<x->R[i]) mx=x->R[i];
  return mx;
}

int Sign(double num)
{
	if(num>0) return 1;
	else if(num<0) return -1;
	else return 0;
}

void Sum(Matrix *out,Matrix *mat)
{
  int i,j,brow,pos;

  brow=0;
  for(i=0;i<mat->m;i++)
  {
    pos=brow;out->R[i]=0;
    for(j=0;j<mat->n;j++)
	{
		out->R[i]+=mat->R[pos];
		pos+=mat->m;
	}
    brow+=1;
  }
	 
}

/*------------ Newton Method ----------------------------------------------*/

struct NewtonStruc
{
  int n;
  Matrix Jac,fx,rhs,x_new,px1,px2;
}
  NewtonData;

void Newton_init(int n)
{
	NewtonData.n=n;
    AllocMatrix(&NewtonData.Jac,n,n);
	AllocMatrix(&NewtonData.fx,n,1);
	AllocMatrix(&NewtonData.x_new,n,1);
	AllocMatrix(&NewtonData.rhs,n,1);
    AllocMatrix(&NewtonData.px1,n,1);
	AllocMatrix(&NewtonData.px2,n,1);
}

void Jacobian(Matrix *Jac,UserFunc Func, double *FuncData, 
			  Matrix *x, Matrix *px1, Matrix *px2, Matrix *DeriSteps)
{
	int n,j,k,ind;
    double x0,dx;

	n=Jac->n;CopyMat(px1,x);CopyMat(px2,x);ind=0;

    for(j=0;j<n;j++)
	{
		x0=x->R[j];dx=DeriSteps->R[j];
        px1->R[j]=x0-dx;
        px2->R[j]=x0+dx;

		for(k=1;k<=n;k++)
		{
            Jac->R[ind]=((*Func)(px2,k,FuncData)-
				         (*Func)(px1,k,FuncData))/(2*dx);        
		    ind++;
		}

        px1->R[j]=x0;
        px2->R[j]=x0;
	}
 
}

int Newton(Matrix *x,
           UserFunc Func, double *FuncData,Matrix *x0,Matrix *DeriSteps,
		   int JacStep,int MaxTime,double Tol)

/* 
   Classical Newton method
 
   DeriSteps - derivative steps (different for each component of x)
   JacStep - re-compute Jacobian only after JacStep steps
   MaxTime - number of maximum allowable iterations
   Tol - x is the solution if ||Func(x,FuncData)||_\infty<=Tol
 
   Newton return 1 if the solution x has been found  */

{
   int i,iter,njac,n;

   n=NewtonData.n;iter=0;njac=JacStep;

   for(i=0;i<n;i++) {NewtonData.fx.R[i]=(*Func)(x0,i+1,FuncData);}
   
   CopyMat(x,x0);

   while(!(InfNorm(&NewtonData.fx)<=Tol) && iter<MaxTime)
   {
	  
      if(njac==JacStep)
	  {
          njac=1;

          Jacobian(&NewtonData.Jac,Func,FuncData, 
			       x,&NewtonData.px1,&NewtonData.px2,DeriSteps);

      } else njac++;

        MatByVec(&NewtonData.rhs,&NewtonData.Jac,x);
        SubMat(&NewtonData.rhs,&NewtonData.fx);

		if(SolveLinearSystem(&NewtonData.x_new,&NewtonData.Jac,&NewtonData.rhs))
		   CopyMat(x,&NewtonData.x_new);
		else
            break;

        for(i=0;i<n;i++) NewtonData.fx.R[i]=(*Func)(x,i+1,FuncData);
		iter++;

   }

   return InfNorm(&NewtonData.fx)<=Tol;
}

void Newton_destroy()
{
   DestroyMatrix(&NewtonData.Jac);
   DestroyMatrix(&NewtonData.fx);
   DestroyMatrix(&NewtonData.x_new);
   DestroyMatrix(&NewtonData.rhs);
   DestroyMatrix(&NewtonData.px1);
   DestroyMatrix(&NewtonData.px2);
}

/*------------------------------ End Newton --------------------------------*/

/*------------------------------- Main Alg ---------------------------------*/


void GetCubeVertices(BoolMatrix *Vertices,unsigned long NVertices,int n)
{
   unsigned long bdigit,mask,bit_on;
   unsigned long ind,i;
   int j;

   ind=0;bdigit=0;
   for(i=0;i<NVertices;i++)
   {
	   for(j=0;j<n;j++)
	   {   
		   mask=(long)pow(2,(double)j);
           bit_on=bdigit & mask;
           if(bit_on) Vertices->R[ind++]=1; else Vertices->R[ind++]=0;
	   }

	   bdigit++;
   }

}


void MainAlg_init()
{
	 int n,i,j,ind;
	 double dx,axis;

	 n=Input.SpaceDim;

	 Work.NVertices=(int)pow(2,(double)n);
	 Work.NCubes=(int)pow((double)Input.Density,(double)n);

     AllocBoolMatrix(&Work.Vertices,n,Work.NVertices);
	 AllocMatrix(&Work.XGrid,Input.Density+1,n);
     AllocBoolMatrix(&Work.Cube,n,1);
	 AllocMatrix(&Work.CubeVert,n,Work.NVertices);
	 AllocMatrix(&Work.Range,Work.NVertices,1);
	 AllocMatrix(&Work.CubeCenter,n,1);
	 AllocMatrix(&Work.Sol,n,Input.SolSize);
	 AllocMatrix(&Work.x,n,1);

	 /* Precompute n-cube vertices */
     GetCubeVertices(&Work.Vertices,Work.NVertices,n);
     
	 ind=0;
	 for(i=0;i<n;i++)
	 {
	   dx=(Input.Limits.R[i+Input.Limits.m]-Input.Limits.R[i])/(double)Input.Density;
	   axis=Input.Limits.R[i]+dx;

       Work.XGrid.R[ind++]=Input.Limits.R[i];
       for(j=1;j<Input.Density;j++) {Work.XGrid.R[ind++]=axis;axis+=dx;}
       Work.XGrid.R[ind++]=Input.Limits.R[i+Input.Limits.m];
	 }

}

void GridAnalysis()
{
   long count;
   int i,j,ind,n,col,MaxGridMark,Outcast,
	   PCubes,NFailed,RSol,SNum,Flag;
   double VMult;

   PCubes=0;NFailed=0;SNum=0;RSol=0;

   VMult=1/(double)Work.NVertices;

   n=Input.SpaceDim;MaxGridMark=Input.Density-1;
   for(i=0;i<n;i++) Work.Cube.R[i]=0;

   for(count=0;count<Work.NCubes;count++)  /* Main cycle - pass through all cubes */
   {
	  /* Translate precomputed abstract cube vertices into real ones */
	  ind=0;
      for(i=0;i<Work.NVertices;i++)
	  {
		 for(j=0;j<n;j++)
		 {
		  col=0;
	      if(Work.Vertices.R[ind]) 
			 Work.CubeVert.R[ind++]=Work.XGrid.R[col+Work.Cube.R[j]+1];
          else
		     Work.CubeVert.R[ind++]=Work.XGrid.R[col+Work.Cube.R[j]];

		  col+=Work.XGrid.m;
		 }
	  }

      /* Outcasting cubes  that cannot contain intersection
         of all nonlinear manifolds */

      Outcast=0;

	  for(i=0;i<n;i++)
	  {
		Work.CubeVert.dis=0;
        for(j=0;j<Work.NVertices;j++) 
		{	
		  Work.Range.R[j]=(*Input.Func)(&Work.CubeVert,i+1,Input.FuncData);
		  Work.CubeVert.dis+=Work.CubeVert.m;
		}
        if(Sign(Min(&Work.Range))==Sign(Max(&Work.Range))) {Outcast=1;break;}
      }

      if(!Outcast)
	  {
         PCubes++;

		 Sum(&Work.CubeCenter,&Work.CubeVert);
		 MatByNum(&Work.CubeCenter,VMult);

		 if(Newton(&Work.x,Input.Func,Input.FuncData,&Work.CubeCenter,&Input.DeriSteps,
		           Input.JacStep,Input.MaxTime,Input.NewtonZero))
		 {
            Work.Sol.dis=0;Flag=1;
            for(i=0;i<SNum;i++) 
			if(Contain(&Work.Sol,&Work.x,Input.SamePoints)) 
			{
				RSol++;Flag=0;break;
			}
			else Work.Sol.dis+=Work.Sol.m;

            if(Flag) {SNum++;AddColumn(&Work.Sol,&Work.x);}
			
		 }
		 else NFailed++;

	  }

      /* Moving to the next cube */
	  for(i=0;i<n;i++)
	  if(Work.Cube.R[i]<MaxGridMark) {Work.Cube.R[i]+=1;break;}
	  else Work.Cube.R[i]=0;

   } /* End main cycle */

   if(SNum) CopyMat(&Output.Sol,&Work.Sol);

   Output.SolutionsNumber.R[0]=SNum;
   Output.ProcessedCubes.R[0]=(double)PCubes;
   Output.AllCubes.R[0]=(double)Work.NCubes;
   Output.NewtonFailed.R[0]=(double)NFailed;
   Output.RedundantSolutions.R[0]=(double)RSol;
}

void MainAlg_destroy()
{
  	 DestroyBoolMatrix(&Work.Vertices);
	 DestroyMatrix(&Work.XGrid);
     DestroyBoolMatrix(&Work.Cube);
     DestroyMatrix(&Work.CubeVert);
	 DestroyMatrix(&Work.Range);
	 DestroyMatrix(&Work.CubeCenter);
	 DestroyMatrix(&Work.Sol);
	 DestroyMatrix(&Work.x);
}

/*----------------------------- End Main Alg ---------------------------------*/

void fill_input_data(mxArray *ptr)
{

     Field2Matrix(ptr,&Input.Limits,"Limits");
	 Field2Matrix(ptr,&Input.DeriSteps,"DeriSteps");
     Input.MaxTime=(int)Field2Double(ptr,"MaxTime");
	 Input.NewtonZero=Field2Double(ptr,"NewtonZero");
	 Input.JacStep=(int)Field2Double(ptr,"JacStep");
	 Input.Density=(int)Field2Double(ptr,"Density");
	 Input.SolSize=(int)Field2Double(ptr,"SolSize");
	 Input.SamePoints=Field2Double(ptr,"SamePoints");
	 Input.Radius=Field2Double(ptr,"Radius");
	 Input.SpaceDim=Input.Limits.m;

	 Input.Func=&kubichek;
	 Input.FuncData=&Input.Radius;
}


void create_work_data()
{
	 dgesv_init(Input.SpaceDim);
	 Newton_init(Input.SpaceDim);
	 MainAlg_init();
}


void create_out_data(mxArray *plhs[])
{
	 mxArray *ptr;
	 char *OutNames[]={"Sol","SolutionsNumber",
		               "AllCubes","ProcessedCubes","NewtonFailed",
					   "RedundantSolutions"};

	 ptr=mxCreateStructMatrix(1,1,6,OutNames);
     CreateOutField(ptr,&Output.Sol,Input.SpaceDim,Input.SolSize,OutNames[0]);
	 CreateOutField(ptr,&Output.SolutionsNumber,1,1,OutNames[1]);
     CreateOutField(ptr,&Output.AllCubes,1,1,OutNames[2]);
     CreateOutField(ptr,&Output.ProcessedCubes,1,1,OutNames[3]);
	 CreateOutField(ptr,&Output.NewtonFailed,1,1,OutNames[4]);
	 CreateOutField(ptr,&Output.RedundantSolutions,1,1,OutNames[5]);

	 plhs[0]=ptr;
}

void do_work()
{
 
     GridAnalysis();     
	 
}

void delete_work_data()
{
	dgesv_destroy();
	Newton_destroy();
	MainAlg_destroy();
}



void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

  if(nrhs<1 || nlhs<1) 
  {
    return;
  } 

  fill_input_data(prhs[0]);
 
  create_work_data();
 
  create_out_data(plhs);
 
  do_work();
 
  delete_work_data();

}
