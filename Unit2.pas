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

  Form2.Chart1.BottomAxis.SetMinMax(0,400+1);//设定坐标
  Form2.Chart1.LeftAxis.SetMinMax(0,36+1);
  Form2.Chart1.Legend.Visible:=true;
  Form2.Chart1.LeftAxis.Labels:=True;
  Form2.Chart1.LeftAxis.AdjustMaxMin;
  Form2.Chart1.LeftAxis.LabelStyle:=talValue;//纵坐标label类型
  Form2.Chart1.LeftAxis.Increment:=1;//设置左坐标的步长为1
  Form2.Chart1.BottomAxis.Increment:=10;
  Form2.Chart1.RemoveAllSeries; //移除chart1的所有series并非删除

  for i:=1 to job_num do
    Tc[i]:=RGB(random(255),random(255),random(255));//工件对应的颜色数组(随机产生)



  for i:=1 to max_i do
    begin

      Ts[i]:=TGanttSeries.Create(Self); //创建Tgantt-series
      Chart1.AddSeries(Ts[i]);//把创建的Tgantt-series放到chart1
      Ts[i].Marks.Visible:=True;//显示每个框图的mark-label
      Ts[i].Marks.Style:=smsLabel;
      Ts[i].Marks.Font.Size:=13;
      Ts[i].Marks.BackColor:=clWhite;//变更marks背景颜色；
      Ts[i].Marks.Font.Style:=[fsItalic]; //设置字体：可以加粗，斜体，下划线~
      Ts[i].Marks.Font.Name:='楷体';
      Ts[i].XValues.DateTime:=False;//禁用X轴的时间显示
      Ts[i].Pointer.HorizSize:=15; //设置甘特图长度
      Ts[i].Pointer.VertSize:=20; //设置甘特图宽度
      Ts[i].Pointer.Style:=psCircle;//甘特图bar的形状


      //Chart1.SeriesList[i-1].Title:='J'+IntToStr(code_arr[i-1])+'-'+'M'+IntToStr(according_machine_order[i-1]);

      {if (process_arr2[according_machine_order[i-1]]=1) then
        begin
          Ts[i].AddGantt(-20.5,0,according_machine_order[i-1],'M'+IntToStr(according_machine_order[i-1]));
        end;}


                 
    end;

end;

end.
