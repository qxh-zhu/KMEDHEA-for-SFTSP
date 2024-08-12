unit Unit2;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls, TeeProcs, TeEngine, Chart , Codec_and,
  Math, Series, GanttCh;

const
  job_num = 60 ;
  max_i = 123 ;

type
  TForm2 = class(TForm)
    btn1: TButton;
    chart1: TChart;
    procedure btn1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form2: TForm2;

implementation

uses Unit1;

{$R *.dfm}

procedure TForm2.btn1Click(Sender: TObject);

var
  i : Integer ;
  Ts : Array [ 1 .. max_i ] of TGanttSeries ;
  Tc : array [ 1 .. Job_num ] of TColor ;

begin

  Form2.Chart1.BottomAxis.SetMinMax(0,400+1);//�趨����
  Form2.Chart1.LeftAxis.SetMinMax(0,36+1);
  Form2.Chart1.Legend.Visible:=true;
  Form2.Chart1.LeftAxis.Labels:=True;
  Form2.Chart1.LeftAxis.AdjustMaxMin;
  Form2.Chart1.LeftAxis.LabelStyle:=talValue;//������label����
  Form2.Chart1.LeftAxis.Increment:=1;//����������Ĳ���Ϊ1
  Form2.Chart1.BottomAxis.Increment:=10;
  Form2.Chart1.RemoveAllSeries; //�Ƴ�chart1������series����ɾ��

  for i:=1 to job_num do
    Tc[i]:=RGB(random(255),random(255),random(255));//������Ӧ����ɫ����(�������)



  for i:=1 to max_i do
    begin

      Ts[i]:=TGanttSeries.Create(Self); //����Tgantt-series
      Chart1.AddSeries(Ts[i]);//�Ѵ�����Tgantt-series�ŵ�chart1
      Ts[i].Marks.Visible:=True;//��ʾÿ����ͼ��mark-label
      Ts[i].Marks.Style:=smsLabel;
      Ts[i].Marks.Font.Size:=13;
      Ts[i].Marks.BackColor:=clWhite;//���marks������ɫ��
      Ts[i].Marks.Font.Style:=[fsItalic]; //�������壺���ԼӴ֣�б�壬�»���~
      Ts[i].Marks.Font.Name:='����';
      Ts[i].XValues.DateTime:=False;//����X���ʱ����ʾ
      Ts[i].Pointer.HorizSize:=15; //���ø���ͼ����
      Ts[i].Pointer.VertSize:=20; //���ø���ͼ���
      Ts[i].Pointer.Style:=psCircle;//����ͼbar����״


      //Chart1.SeriesList[i-1].Title:='J'+IntToStr(code_arr[i-1])+'-'+'M'+IntToStr(according_machine_order[i-1]);

      {if (process_arr2[according_machine_order[i-1]]=1) then
        begin
          Ts[i].AddGantt(-20.5,0,according_machine_order[i-1],'M'+IntToStr(according_machine_order[i-1]));
        end;}


                 
    end;

end;

end.
