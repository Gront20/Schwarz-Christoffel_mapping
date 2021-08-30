//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
extern "C"
{
 void _stdcall FR1(int*,int*);
 int _stdcall FF1(int*);
};
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TButton *Button1;
        TImage *Image1;
        TImage *Image3;
        TImage *Image5;
        TPanel *Panel1;
        TImage *Image4;
        TLabel *Label1;
        TListBox *ListBox1;
        TLabel *Label2;
        TListBox *ListBox2;
        void __fastcall Button1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
 