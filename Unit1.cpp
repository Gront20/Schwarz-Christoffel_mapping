//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
 
#include <math.h>
#include "Unit1.h"

#define COUNT_PER_SHAPE 100

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

typedef void __stdcall SchwarzChristoff_Acessor(int *,double *,double *,double *,int *,double *,double *,double *,double *,double *);
typedef void __stdcall PolygonToDisk(int *,double *,double *,double *,double *,double *,int *,double *,double *,int*,double *);
typedef void __stdcall DiskToPolygon(int *,double *,double *,double *,double *,double *,int *,double *,double *,int*,double *);


//---------------------------------------------------------------------------
struct Complex {
  double re;
  double im;
  Complex(void)
  {
    re=0;
    im=0;
  };
  Complex(double are, double aim)
  {
   re=are;
   im=aim;
  };
  Complex operator -(const Complex& v) const
  {
    return Complex(re-v.re,im-v.im);
  }
  bool operator <(const Complex &v) const
  {
    return re<v.re || (re==v.re && im<v.im);
  }
  double Complex::Abs() const
  {
   return ::sqrt(re * re + im * im);
  }
};

//---------------------------------------------------------------------------
class TStr_SchwarzChristoff
{
public:
  TList *AVal;
  int n;
  double *w;
  double *wc;
  double *betam;
  int nptsq;
  double *tol;
  double *errest;
  double *c;
  double *z;
  double *qwork;
  __fastcall TStr_SchwarzChristoff(TList *FVal);
};

//---------------------------------------------------------------------------
class TStr_PolygonToDisk
{
public:
  int n;
  double *c;
  double *z;
  double *wc;
  double *w;
  double *betam;
  int nptsq;
  double *qwork;
  double *ww;
  int npred;
  double *zz;
  __fastcall TStr_PolygonToDisk(TStr_SchwarzChristoff &val);
};

//---------------------------------------------------------------------------
class TStr_DiskToPolygon
{
public:
  int n;
  double *c;
  double *z;
  double *wc;
  double *w;
  double *betam;
  int nptsq;
  double *qwork;
  double *zz;
  int npred;
  double *ww;
  __fastcall TStr_DiskToPolygon(TStr_SchwarzChristoff &val);
};

//---------------------------------------------------------------------------

class TPolygon
{
private:
  double VneshUgol(Complex *A,Complex *B,Complex *C);
public:
  TList *AVal;
  double *beta;
  __fastcall TPolygon(TList *FVal);
};
//---------------------------------------------------------------------------

__fastcall TStr_SchwarzChristoff::TStr_SchwarzChristoff(TList *FVal)
{
  AVal=new TList();
  for(int i=0;i<FVal->Count;i++)
  {
     Complex *t2=new Complex();
     t2->re=((Complex *)(FVal->Items[i]))->re;
     t2->im=((Complex *)(FVal->Items[i]))->im;
     AVal->Add(t2);
  }
  TPolygon pg(FVal);
  w=new double[2*FVal->Count];
  wc=new double[2];
  betam=new double[FVal->Count];
  nptsq=12;
  tol=new double;
  errest=new double;
  c=new double[2];
  z=new double[2*FVal->Count];
  qwork=new double[nptsq*(2*FVal->Count+3)];
  for(int i=0;i<FVal->Count;i++)
  {
    w[2*i]=((Complex *)(FVal->Items[i]))->re;
    w[2*i+1]=((Complex *)(FVal->Items[i]))->im;
  }
  n=FVal->Count;
  double q1=0,q2=0;
  for(int i=0;i<FVal->Count;i++)
  {
     q1+=((Complex *)(FVal->Items[i]))->re;
     q2+=((Complex *)(FVal->Items[i]))->im;
  }
  wc[0]=q1/FVal->Count;
  wc[1]=q2/FVal->Count;
  for(int i=0;i<FVal->Count;i++)
  {
      betam[i]=-pg.beta[i]/3.141592653589793238462643383279 ;
  }
}

//---------------------------------------------------------------------------

__fastcall TStr_PolygonToDisk::TStr_PolygonToDisk(TStr_SchwarzChristoff &val)
{
  n=val.n;
  c=val.c;
  z=val.z;
  wc=val.wc;
  w=val.w;
  betam=val.betam;
  nptsq=val.nptsq;
  qwork=val.qwork;
}

//---------------------------------------------------------------------------

__fastcall TStr_DiskToPolygon::TStr_DiskToPolygon(TStr_SchwarzChristoff &val)
{
  n=val.n;
  c=val.c;
  z=val.z;
  wc=val.wc;
  w=val.w;
  betam=val.betam;
  nptsq=val.nptsq;
  qwork=val.qwork;
}

//---------------------------------------------------------------------------

__fastcall TPolygon::TPolygon(TList *FVal)
{
  AVal=new TList();
  for(int i=0;i<FVal->Count;i++)
  {
     Complex *t2=new Complex();
     t2->re=((Complex *)(FVal->Items[i]))->re;
     t2->im=((Complex *)(FVal->Items[i]))->im;
     AVal->Add(t2);
  }
  int n = AVal->Count;

  beta=new double[n];
  beta[0] = VneshUgol((Complex *)(AVal->Items[n-1]),(Complex *)(AVal->Items[0]), (Complex *)(AVal->Items[1]));
  beta[n - 1] = VneshUgol((Complex *)(AVal->Items[n-2]),(Complex *)(AVal->Items[n-1]),(Complex *)(AVal->Items[0]));
  for (int i = 1; i < n - 1; ++i) {
    beta[i] = VneshUgol((Complex *)(AVal->Items[i-1]),(Complex *)(AVal->Items[i]),(Complex *)(AVal->Items[i+1]));
  }
}
//---------------------------------------------------------------------------
double TPolygon::VneshUgol(Complex *A,Complex *B,Complex *C) {
 Complex p1=Complex((*B)-(*A));
 Complex p2=Complex((*C)-(*B));
 double t1=p1.re*p2.im-p1.im*p2.re;
 double t2=p1.re*p2.re+p1.im*p2.im;
 return atan2(t1,t2);
}

//---------------------------------------------------------------------------
double __fastcall MaxValuePoint(TList *points)
{
   double m=0;
   for(int i=0;i<points->Count;i++)
   {
      if(fabs(((Complex *)(points->Items[i]))->re)>m) m=fabs(((Complex *)(points->Items[i]))->re);
      if(fabs(((Complex *)(points->Items[i]))->im)>m) m=fabs(((Complex *)(points->Items[i]))->im);
   }
   return m;
}
//---------------------------------------------------------------------------
bool new_inside(double x,double y,TList *polygon)
{
    static const int q_patt[2][2]= { {0,1}, {3,2} };
    if (polygon->Count<3) return false;

    Complex pred_pt=*(Complex *)(polygon->Items[0]);
    pred_pt.re-=x;
    pred_pt.im-=y;

    int pred_q=q_patt[pred_pt.im<0][pred_pt.re<0];

    int w=0;

    for(int i=0;i<polygon->Count;i++)
    {
       Complex *temp=(Complex *)(polygon->Items[i]);
       Complex cur_pt=*temp;
       cur_pt.re-=x;
       cur_pt.im-=y;

       int q=q_patt[cur_pt.im<0][cur_pt.re<0];
       switch (q-pred_q)
       {
           case -3:++w;break;
           case 3:--w;break;
           case -2:if(pred_pt.re*cur_pt.im>=pred_pt.im*cur_pt.re) ++w;break;
           case 2:if(!(pred_pt.re*cur_pt.im>=pred_pt.im*cur_pt.re)) --w;break;
       }
       pred_pt = *temp;
       pred_q = q;
    }
    return w!=0;
}

//---------------------------------------------------------------------------
int inside(double px,double py,TList *lst) {
    int i,j,s;
    int n=lst->Count;
    s=0; j=n-1;
    for(i=0;i<n;j=i++) {
        Complex *t1=(Complex *)(lst->Items[i]);
        double xi=t1->re;
        double yi=t1->im;
        t1=(Complex *)(lst->Items[j]);
        double xj=t1->re;
        double yj=t1->im;
        if ( (py<yi ^ py<yj) || py==yi || py==yj ) {
            if ( (px<xi ^ px<xj) || px==xi || px==xj ) {
                if ( yi<yj ^ (xi-px)*(yj-py)>(yi-py)*(xj-px) ) s^=1;
            } else {
                if (px>xi && yi!=yj) s^=1;
            }
        }
    }
    if(s) s=new_inside(px,py,lst);
    return s;
}

//---------------------------------------------------------------------------
double __fastcall MaxLine(TList *points)
{
   double m=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      int l=i+1;
      if(l>d-1) l=d-i-1;
      Complex *t2=(Complex *)(points->Items[l]);
      Complex t3=*t1-*t2;
      if(m<t3.Abs()) m=t3.Abs();
   }
   return m;
}

//---------------------------------------------------------------------------
void __fastcall Zapol(TList *A,TList *points)
{
    double pntcount=MaxLine(A)/COUNT_PER_SHAPE;
    int L=A->Count;
    for(int i=0;i<L;i++)
    {
      Complex *t1=(Complex *)(A->Items[i]);
      int l=i+1;
      if(l>L-1) l=L-i-1;
      Complex *t2=(Complex *)(A->Items[l]);
      Complex t3=*t1-*t2;
      int n=ceil(t3.Abs()/pntcount);
      for(int k=0;k<n;k++)
      {
        double x=t1->re+k*pntcount*(t2->re-t1->re)/t3.Abs();
        double y=t1->im+k*pntcount*(t2->im-t1->im)/t3.Abs();
        Complex *tmp=new Complex(x,y);
        points->Add(tmp);
      }
    }
}


//---------------------------------------------------------------------------
void DrawShape(TList *vx,TImage *Image_drw,int isline,double pntcount1)
{
//  double pntcount1=MaxValuePoint(vx);
  int mash;
  mash=int(Image_drw->Height/(pntcount1*2));
  int offsetX=Image_drw->Width/2;
  int offsetY=Image_drw->Height/2;
  Image_drw->Canvas->MoveTo(offsetX,0);
  Image_drw->Canvas->LineTo(offsetX,Image_drw->Height);
  Image_drw->Canvas->MoveTo(0,offsetY);
  Image_drw->Canvas->LineTo(Image_drw->Width,offsetY);
  int l;
  l=vx->Count;
  for(int i=0;i<l;i++)
  {
         int x1=int(((Complex *)(vx->Items[i]))->re*mash+double(offsetX));
         int y1=int(double(offsetY)-((Complex *)(vx->Items[i]))->im*mash);
         if(!isline) Image_drw->Canvas->Ellipse(x1,y1,x1+2,y1+2);
         else
         {
            if(i<l-1)
            {
              int x2=int(((Complex *)(vx->Items[i+1]))->re*mash+double(offsetX));
              int y2=int(double(offsetY)-((Complex *)(vx->Items[i+1]))->im*mash);
              Image_drw->Canvas->MoveTo(x1,y1);
              Image_drw->Canvas->LineTo(x2,y2);
            } else
            {
              int x2=int(((Complex *)(vx->Items[0]))->re*mash+double(offsetX));
              int y2=int(double(offsetY)-((Complex *)(vx->Items[0]))->im*mash);
              Image_drw->Canvas->MoveTo(x1,y1);
              Image_drw->Canvas->LineTo(x2,y2);
            }
         }
         Application->ProcessMessages();
  }
}
//---------------------------------------------------------------------------
double __fastcall Razmer(TList *points)
{
   double m=0;
   int d=points->Count;
   for(int i=0;i<d;i++)
   {
      Complex *t1=(Complex *)(points->Items[i]);
      for(int j=0;j<d;j++)
      {
         Complex *t2=(Complex *)(points->Items[j]);
         double v=sqrt(pow((t1->re-t2->re),2)+pow((t1->im-t2->im),2));
         if(m<v) m=v;
      }
   }
   return m;
}

//---------------------------------------------------------------------------

void GenerateGrid(TList *FVal,TList *res)
{
     TList *_z=new TList();
     for(int i=0;i<res->Count;i++) delete res->Items[i];
     res->Clear();

     Zapol(FVal,_z);
     double pntcount=Razmer(_z);
     double xstep=2*pntcount/100;
     double ystep=2*pntcount/100;

     double hx=2*pntcount/400;
     double hy=2*pntcount/400;

     for(int i=0;i<100;i++)
     {
        double x=pntcount-xstep*i;
        for(int j=0;j<400;j++)
        {
           double y=pntcount-hy*j;
           Complex *pt=new Complex(x,y);
           if(inside(pt->re,pt->im,_z)) res->Add(pt);
        }
     }
     for(int i=0;i<100;i++)
     {
        double y=pntcount-ystep*i;
        for(int j=0;j<400;j++)
        {
           double x=pntcount-hx*j;
           Complex *pt=new Complex(x,y);
           if(inside(pt->re,pt->im,_z)) res->Add(pt);
        }
     }
     DrawShape(_z,Form1->Image1,0,pntcount);
     DrawShape(res,Form1->Image1,0,pntcount);     
     for(int i=0;i<_z->Count;i++) delete _z->Items[i];
     _z->Clear();
     delete _z;
}



//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
  HMODULE  hDll ;

  Button1->Enabled=false;
  hDll = LoadLibrary("sclib.dll") ;
  if ( hDll == NULL )
  {
      ShowMessage("Íå âîçìîæíî çàãðóçèòü dll ! "+IntToStr(GetLastError()));
      return;
  }
  SchwarzChristoff_Acessor *pSchwarzChristoff = (SchwarzChristoff_Acessor *)GetProcAddress(hDll,"schwchris_" ) ;
  if(pSchwarzChristoff==NULL )
  {
      ShowMessage("Àäðåñ ïðîöåäóðû SchwChris íå íàéäåí");
      return ;
  }
  PolygonToDisk *pPolygonToDisk = (PolygonToDisk *)GetProcAddress(hDll,"poltodisk_" ) ;
  if(pPolygonToDisk==NULL )
  {
      ShowMessage("Àäðåñ ïðîöåäóðû polygtodisk íå íàéäåí");
      return ;
  }

  DiskToPolygon *pDiskToPolygon = (DiskToPolygon *)GetProcAddress(hDll,"disktopol_" ) ;
  if(pDiskToPolygon==NULL )
  {
      ShowMessage("Àäðåñ ïðîöåäóðû disktopolyg íå íàéäåí");
      return ;
  }


//===================================================================================

  TList *vx=new TList;

  //Äèíîçàâð
  Complex *t1=new Complex(3,-1);
  t1=new Complex(3,-1);
  vx->Add(t1);
  t1=new Complex(1,2);
  vx->Add(t1);
  t1=new Complex(1,4);
  vx->Add(t1);
  t1=new Complex(3,2);
  vx->Add(t1);
  t1=new Complex(4,3);
  vx->Add(t1);
  t1=new Complex(2,6);
  vx->Add(t1);
  t1=new Complex(-1,6);
  vx->Add(t1);
  t1=new Complex(-2,4);
  vx->Add(t1);
  t1=new Complex(-2,2);
  vx->Add(t1);
  t1=new Complex(-3,3);
  vx->Add(t1);
  t1=new Complex(-3,0);
  vx->Add(t1);
  t1=new Complex(-2,1);
  vx->Add(t1);
  t1=new Complex(-2,-3);
  vx->Add(t1);
  t1=new Complex(-4,-2);
  vx->Add(t1);
  t1=new Complex(-2,-4);
  vx->Add(t1);
  t1=new Complex(1,-4);
  vx->Add(t1);

  TStr_SchwarzChristoff ssm(vx);
//  DrawShape(vx,Image1,1);

//Âû÷èñëåíèå àöåññîðíûõ ïàðàìåòðîâ îòîáðàæåíèÿ äèíîçàâðà

//===================================================================================
  (*pSchwarzChristoff)(&(ssm.n),ssm.w,ssm.wc,ssm.betam,&(ssm.nptsq),ssm.tol,ssm.errest,ssm.c,ssm.z,ssm.qwork);
//===================================================================================
  ListBox1->Items->Add("wc="+FloatToStrF(ssm.wc[0],ffFixed,7,4)+" + i"+FloatToStrF(ssm.wc[1],ffFixed,7,4));
  ListBox1->Items->Add("C="+FloatToStrF(ssm.c[0],ffFixed,7,4)+" + i"+FloatToStrF(ssm.c[1],ffFixed,7,4));
  for(int i=0;i<vx->Count;i++)
  {
     if(ssm.z[2*i+1]>=0)  ListBox1->Items->Add("z["+IntToStr(i)+"]="+FloatToStrF(ssm.z[2*i],ffFixed,7,4)+" + i"+FloatToStrF(ssm.z[2*i+1],ffFixed,7,4));
     else ListBox1->Items->Add("z["+IntToStr(i)+"]="+FloatToStrF(ssm.z[2*i],ffFixed,7,4)+" - i"+FloatToStrF(fabs(ssm.z[2*i+1]),ffFixed,7,4));
  }


  TStr_PolygonToDisk spd(ssm);
  TList *vv=new TList();
  GenerateGrid(vx,vv);

  TStringList *tst=new TStringList();
  for(int i=0;i<vx->Count;i++)
  {
    tst->Add("z["+IntToStr(i)+"]=x"+FloatToStr(ssm.z[2*i])+"; y="+FloatToStr(ssm.z[2*i+1]));
  }
  tst->SaveToFile("out_z.txt");
  delete tst;


  tst=new TStringList();
  for(int i=0;i<vv->Count;i++)
  {
    tst->Add("i="+IntToStr(i)+" x="+FloatToStr(((Complex *)(vv->Items[i]))->re)+"; y="+FloatToStr(((Complex *)(vv->Items[i]))->im));
  }
  tst->SaveToFile("out1.txt");
  delete tst;

//  double pntcount=MaxValuePoint(vv);
//  DrawShape(vv,Image1,0,pntcount);
  spd.ww=new double[2*vv->Count];
  for(int i=0;i<vv->Count;i++)
  {
    spd.ww[2*i]=((Complex *)(vv->Items[i]))->re;
    spd.ww[2*i+1]=((Complex *)(vv->Items[i]))->im;
  }
  spd.npred=vv->Count;
  spd.zz=new double[2*spd.npred];


//Êîíôîðìíîå îòîáðàæåíèå âíóòðåííîñòè äèíîçàâðà íà åäèíè÷íûé äèñê

//===================================================================================
    (*pPolygonToDisk)(&(spd.n),spd.c,spd.z,spd.wc,spd.w,spd.betam,&(spd.nptsq),spd.qwork,spd.ww,&(spd.npred),spd.zz);
//===================================================================================

  TList *vr=new TList();
  for(int i=0;i<spd.npred;i++)
  {
      if(::sqrt(spd.zz[2*i]*spd.zz[2*i]+spd.zz[2*i+1]*spd.zz[2*i+1])<1)
      {
        Complex *t1=new Complex(spd.zz[2*i],spd.zz[2*i+1]);
        vr->Add(t1);
      }
  }

  tst=new TStringList();
  for(int i=0;i<vr->Count;i++)
  {
    tst->Add("i="+IntToStr(i)+" x="+FloatToStr(((Complex *)(vr->Items[i]))->re)+"; y="+FloatToStr(((Complex *)(vr->Items[i]))->im));
  }
  tst->SaveToFile("out2.txt");
  delete tst;
  double pntcount1=MaxValuePoint(vr);
  DrawShape(vr,Image5,0,pntcount1);

  TList *vx1=new TList();


//Êîò
  t1=new Complex(3,-4);
  vx1->Add(t1);
  t1=new Complex(4,0);
  vx1->Add(t1);
  t1=new Complex(5,1);
  vx1->Add(t1);
  t1=new Complex(5,2);
  vx1->Add(t1);
  t1=new Complex(4,3);
  vx1->Add(t1);
  t1=new Complex(3,2);
  vx1->Add(t1);
  t1=new Complex(2,0);
  vx1->Add(t1);
  t1=new Complex(1,5);
  vx1->Add(t1);
  t1=new Complex(-1,7);
  vx1->Add(t1);
  t1=new Complex(-1,8);
  vx1->Add(t1);
  t1=new Complex(0,9);
  vx1->Add(t1);
  t1=new Complex(0,11);
  vx1->Add(t1);
  t1=new Complex(-1,12);
  vx1->Add(t1);
  t1=new Complex(-2,12);
  vx1->Add(t1);
  t1=new Complex(-3,11);
  vx1->Add(t1);
  t1=new Complex(-3,9);
  vx1->Add(t1);
  t1=new Complex(-2,8);
  vx1->Add(t1);
  t1=new Complex(-2,7);
  vx1->Add(t1);
  t1=new Complex(-6,3);
  vx1->Add(t1);
  t1=new Complex(-6,-4);
  vx1->Add(t1);
  t1=new Complex(-5,-4);
  vx1->Add(t1);
  t1=new Complex(-5,2);
  vx1->Add(t1);
  t1=new Complex(-4,3);
  vx1->Add(t1);
  t1=new Complex(-2,1);
  vx1->Add(t1);
  t1=new Complex(-2,-1);
  vx1->Add(t1);
  t1=new Complex(-3,0);
  vx1->Add(t1);
  t1=new Complex(-4,-1);
  vx1->Add(t1);
  t1=new Complex(-4,-4);
  vx1->Add(t1);

  double pntcount2=MaxValuePoint(vx1);
  DrawShape(vx1,Image3,1,pntcount2);
  TStr_SchwarzChristoff ssm1(vx1);

//Âû÷èñëåíèå àöåññîðíûõ ïàðàìåòðîâ îòîáðàæåíèÿ êîòà

//===================================================================================
  (*pSchwarzChristoff)(&(ssm1.n),ssm1.w,ssm1.wc,ssm1.betam,&(ssm1.nptsq),ssm1.tol,ssm1.errest,ssm1.c,ssm1.z,ssm1.qwork);
//===================================================================================
  ListBox2->Items->Add("wc="+FloatToStrF(ssm1.wc[0],ffFixed,7,4)+" + i"+FloatToStrF(ssm1.wc[1],ffFixed,7,4));
  ListBox2->Items->Add("C="+FloatToStrF(ssm1.c[0],ffFixed,7,4)+" + i"+FloatToStrF(ssm1.c[1],ffFixed,7,4));
  for(int i=0;i<vx1->Count;i++)
  {
     if(ssm1.z[2*i+1]>=0)  ListBox2->Items->Add("z["+IntToStr(i)+"]="+FloatToStrF(ssm1.z[2*i],ffFixed,7,4)+" + i"+FloatToStrF(ssm1.z[2*i+1],ffFixed,7,4));
     else ListBox2->Items->Add("z["+IntToStr(i)+"]="+FloatToStrF(ssm1.z[2*i],ffFixed,7,4)+" - i"+FloatToStrF(fabs(ssm1.z[2*i+1]),ffFixed,7,4));
  }



//Ïîäñòàíîâêà äàííûõ îòîáðàæåíèÿ äèíîçàâðà íà åäèíè÷íûé äèñê
//===================================================================================
  TStr_DiskToPolygon sdp(ssm1);
  sdp.zz=spd.zz;
  sdp.npred=spd.npred;
  sdp.ww=new double[2*spd.npred];
//===================================================================================


//Êîíôîðìíîå îòîáðàæåíèå åäèíè÷íîãî äèñêà, ñîäåðæàùåãî îòîáðàæåíèå âíóòðåííîñòè äèíîçàâðà âî âíóòðåííîñòü îãðàíè÷åííîé ïîëèãîíîì , îïèñûâàþùåãî êîòà

//===================================================================================
  (*pDiskToPolygon)(&(sdp.n),sdp.c,sdp.z,sdp.wc,sdp.w,sdp.betam,&(sdp.nptsq),sdp.qwork,sdp.zz,&(sdp.npred),sdp.ww);
//===================================================================================

  TList *vr1=new TList();
  for(int i=0;i<sdp.npred;i++)
  {
    Complex pt(sdp.ww[2*i],sdp.ww[2*i+1]);
    if(inside(pt.re,pt.im,vx1))
    {
       Complex *t1=new Complex(sdp.ww[2*i],sdp.ww[2*i+1]);
       vr1->Add(t1);
    }
  }
  tst=new TStringList();
  for(int i=0;i<vr1->Count;i++)
  {
    tst->Add("i="+IntToStr(i)+" x="+FloatToStr(((Complex *)(vr1->Items[i]))->re)+"; y="+FloatToStr(((Complex *)(vr1->Items[i]))->im));
  }
  tst->SaveToFile("out3.txt");
  delete tst;

  double pntcount3=MaxValuePoint(vx1);
  int mash;
  mash=int(Image4->Height/(pntcount3*2));
  int offsetX=Image4->Width/2;
  int offsetY=Image4->Height/2;
  Image4->Canvas->MoveTo(offsetX,0);
  Image4->Canvas->LineTo(offsetX,Image4->Height);
  Image4->Canvas->MoveTo(0,offsetY);
  Image4->Canvas->LineTo(Image4->Width,offsetY);
  int l;
  l=vx1->Count;
  for(int i=0;i<l;i++)
  {
         int x1=int(((Complex *)(vx1->Items[i]))->re*mash+double(offsetX));
         int y1=int(double(offsetY)-((Complex *)(vx1->Items[i]))->im*mash);
         if(i<l-1)
         {
              int x2=int(((Complex *)(vx1->Items[i+1]))->re*mash+double(offsetX));
              int y2=int(double(offsetY)-((Complex *)(vx1->Items[i+1]))->im*mash);
              Image4->Canvas->MoveTo(x1,y1);
              Image4->Canvas->LineTo(x2,y2);
          } else
          {
              int x2=int(((Complex *)(vx1->Items[0]))->re*mash+double(offsetX));
              int y2=int(double(offsetY)-((Complex *)(vx1->Items[0]))->im*mash);
              Image4->Canvas->MoveTo(x1,y1);
              Image4->Canvas->LineTo(x2,y2);
          }
          Application->ProcessMessages();
  }
  l=vr1->Count;
  for(int i=0;i<l;i++)
  {
         int x1=int(((Complex *)(vr1->Items[i]))->re*mash+double(offsetX));
         int y1=int(double(offsetY)-((Complex *)(vr1->Items[i]))->im*mash);
         Image4->Canvas->Ellipse(x1,y1,x1+2,y1+2);
         Application->ProcessMessages();
  }
}
//---------------------------------------------------------------------------
