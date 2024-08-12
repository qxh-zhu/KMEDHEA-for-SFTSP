unit Unit1;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Codec_and, math, TeEngine, Series, ExtCtrls, TeeProcs, Chart;

const
  max_i= 123 ;

type
  TForm1 = class(TForm)
    btn1: TButton;
    KMEDHEA: TButton;
    btn3: TButton;
    btn4: TButton;
    HHMEDA_compare: TButton;
    NFOA: TButton;
    HEDA: TButton;
    CEDA: TButton;
    ceda1: TButton;
    KMEA: TButton;
    CCIWO: TButton;
    CSRS: TButton;

    procedure btn1Click(Sender: TObject);
    procedure KMEDHEAClick(Sender: TObject);
    procedure btn3Click(Sender: TObject);
    procedure btn4Click(Sender: TObject);
    procedure HHMEDA_compareClick(Sender: TObject);
    procedure NFOAClick(Sender: TObject);
    procedure HEDAClick(Sender: TObject);
    procedure CEDAClick(Sender: TObject);
    procedure ceda1Click(Sender: TObject);
    procedure KMEAClick(Sender: TObject);
    procedure CCIWOClick(Sender: TObject);
    procedure CSRSClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }

  end;



var
  Form1: TForm1;

    test , hand : arr ;
    accessory : arr ;
    setup_time : two_arr1 ;
    o_arr : arr ;
    M_mt_arr : three_arr ;
    P_mt_arr : three_arr ;
    source_match : two_arr4 ;

implementation

uses Unit2 ;

{$R *.dfm}

procedure TForm1.btn1Click(Sender: TObject);

var
    strL : TStringList ;
    str, str_P : string ;
    n , m , i , Pos_KongGe , str_length , GJ , GX : integer ;


begin

  strL := TStringList.create ;

  strL.LoadFromFile ( './dataset/' + '10' + '.txt' ) ;

  fillchar ( test , sizeof( arr ) , 0 ) ;

  fillchar ( hand , sizeof( arr ) , 0 ) ;

  fillchar ( accessory , sizeof( arr ) , 0 ) ;

  fillchar ( setup_time , sizeof( two_arr1 ) , 0 ) ;

  fillchar ( M_mt_arr , sizeof( three_arr ) , 0 ) ;

  fillchar ( P_mt_arr , sizeof( three_arr ) , 0 ) ;

  for n := 1 to 3 do                                                            //测试资源数目

  begin

    if n = 1 then

    begin

      str := strL [ n - 1 ] ;

      str_length := length(str) ;

      for m := 1 to 3 do

      begin

        Pos_KongGe := pos( ' ' , str ) ;

        str_P := copy( str , 1 , Pos_KongGe-1 ) ;

        test [ m ] := strtoint ( str_P ) ;

        if m <> 3 then

           str := copy ( str , Pos_KongGe + 1 , str_length - Pos_KongGe ) ;

      end ;

    end ;

    if n = 2  then                                                              //处理资源数目

    begin

      str := strL [ n - 1 ] ;

      str_length := length ( str ) ;

      for m := 1 to 3 do

      begin

        Pos_KongGe := pos ( ' ' , str ) ;

        str_P := copy ( str , 1 , Pos_KongGe-1 ) ;

        hand [ m ] := strtoint ( str_P ) ;

        if m <> 3 then

           str := copy ( str , Pos_KongGe + 1 , str_length - Pos_KongGe ) ;

      end ;

    end ;

    if n = 3  then
                                                                                //附件资源数目
    begin

      str := strL [ n - 1 ] ;

      str_length := length ( str ) ;

      for m := 1 to 4 do

      begin

        Pos_KongGe := pos ( ' ' , str ) ;

        str_P := copy ( str , 1 , Pos_KongGe-1 ) ;

        accessory [ m ] := strtoint ( str_P ) ;

        if m <> 4 then

           str := copy ( str , Pos_KongGe + 1 , str_length - Pos_KongGe ) ;

      end ;

    end ;

  end;

  for n := 7 to 42 do                                                           //资源配置

  begin

    str := strL [ n - 1 ] ;

    str_length := length ( str ) ;

    for m := 1 to 3 do

    begin

      Pos_KongGe := pos ( ' ' , str ) ;

      str_P := copy ( str , 1 , Pos_KongGe-1 ) ;

      source_match [ n - 6 , m ] := strtoint ( str_P ) ;

      if m <> 3 then

         str := copy ( str , Pos_KongGe + 1 , str_length - Pos_KongGe ) ;

    end ;

  end ;

  for n := 44 to 79 do
                                                                                //设置时间
  begin

    str := strL [ n - 1 ] ;

    str_length := length ( str ) ;

    for m := 1 to 36 do

    begin

      Pos_KongGe := pos ( ' ' , str ) ;

      str_P := copy ( str , 1 , Pos_KongGe-1 ) ;

      setup_time [ n - 43 , m ] := strtoint ( str_P ) ;

      if m <> 36 then

         str := copy ( str , Pos_KongGe + 1, str_length - Pos_KongGe ) ;

    end ;

  end ;

  i := 1 ;

  for n := 81 to 80 + max_i do                                                  //机器 , 加工时间矩阵

  begin

    str := strL [ n - 1 ] ;

    str_length := length ( str ) ;


    m := 1 ;

    repeat

      Pos_KongGe := pos ( ' ' , str ) ;

      str_P := copy ( str , 1 , Pos_KongGe - 1) ;

      if m = 1 then

      begin

        GJ := strtoint ( str_P ) ;

        o_arr [ i ] := GJ ;

        i := i + 1 ;

      end ;

      if m = 2 then

      begin

        GX := strtoint ( str_P ) ;

      end ;

      if m = 3 then

      begin

        M_mt_arr [ GJ , GX , 1 ] := strtoint ( str_P ) + 1 ;

      end;

      if m = 4 then

      begin

        P_mt_arr [ GJ , GX , 1 ] := strtoint ( str_P ) ;

      end;

      if m = 5 then

      begin

        M_mt_arr [ GJ , GX , 2 ] := strtoint ( str_P ) + 1 ;

      end;

      if m = 6 then

      begin

        P_mt_arr [ GJ , GX , 2 ] := strtoint ( str_P ) ;

      end;

      if m = 7 then

      begin

        M_mt_arr [ GJ , GX , 3 ] := strtoint ( str_P ) + 1 ;

      end;

      if m = 8 then

      begin

        P_mt_arr [ GJ , GX , 3 ] := strtoint ( str_P ) ;

      end;

      str := copy ( str , Pos_KongGe + 1 , str_length - Pos_KongGe ) ;

      m := m + 1 ;

    until str = '' ;

  end ;

end ;


procedure TForm1.KMEDHEAClick(Sender: TObject) ;


begin

  main ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end ;



procedure TForm1.btn3Click(Sender: TObject);
begin

  HHEDA ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.btn4Click(Sender: TObject);
begin

  QHH ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.HHMEDA_compareClick(Sender: TObject);
begin

  HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;
end;

procedure TForm1.NFOAClick(Sender: TObject);
begin

  HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.HEDAClick(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.CEDAClick(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.ceda1Click(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.KMEAClick(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.CCIWOClick(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

procedure TForm1.CSRSClick(Sender: TObject);
begin
    HHMEDA_compare1 ( o_arr ,  test , hand , accessory ,           //工序，三种资源数
         source_match ,                               //所需资源类型
         setup_time ,                                 //设置时间
         M_mt_arr , P_mt_arr  ) ;

end;

end.
