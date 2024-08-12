unit Codec_and;

interface
  uses
    Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
    Dialogs, StdCtrls, math, TeEngine, Series, ExtCtrls, TeeProcs, Chart;

const
     max_i = 123 ;                              //工序数
     initially_generated_individuals = 100 ;    //初始生成的个体数量
     hls_popsize = 30 ;                         //高层操作种群大小
     hls_length = 10 ;                          //高层操作序长度   也是高层种群取优数量
     MC_length = hls_length - 1 ;               //高层序列数-1为矩阵长
     MC_width_and_hight = 8 ;                   //矩阵宽和高
     r = 0.4 ;                                 //学习率

     QHH_OP = 8 ;                               //QHH低层算子数量

type
      arr = array [ 1 .. MAX_i ] of Integer ;   //arr：长度为总工序数的数组
      arrs = array [ 1 .. 100 ] of Integer ;
      arr_loop = array [ 1 .. 36 ] of Integer ; //loop:用来转化前插解码后甘特图
      S_arr = array [ 1 .. 3 ] of Integer ;
      two_arr1 = array [ 1 .. 36 ] of array [ 1 .. 36 ] of Integer ;
      two_arr2 = array [ 1 .. 100 ] of array [ 1 .. 3 ] of Integer ;
      two_arr3 = array [ 1 .. 36 ] of array [ 1 .. 2 ] of Integer ;
      two_arr4 = array [ 1 .. 36 ] of array [ 1 .. 3 ] of Integer ;
      three_arr1 = array [ 1 .. 36 ] of array [ 1 .. 36 ] of array [ 1 .. 2 ] of Integer ;
      three_arr = array [ 1 .. 100 ] of array [ 1 .. 3 ] of array [ 1 .. 3 ] of Integer ;
      psum_arr = array [ 1 .. MC_width_and_hight ] of Real ;    //概率累计数列

type  pro_pos = array [ 1 .. 36 ] of two_arr3 ;

type Indiv1 = record              //Indiv1 ： 工序序 ； 机器序 ； 加工时间序
              op_arr : arr ;
              ma_arr : arr ;
              pr_arr : arr ;
              end  ;

type SBZM = record                                        // SBZM：
           pri_op : array [ 1.. max_i ] of integer ;      { 原工序 }
           op : array [ 1.. max_i ] of integer ;          { 经过规则生成的工序’ }
           ma : array [ 1.. max_i ] of integer ;          { 相匹配的机器序’ }
           fitness : Double ;                             { 完工时间 }
           end ;

type hls_indiv = record                                    //高层个体
           arr : array [ 1.. hls_length ] of integer ;     { 一个个算子 }
           evalue : double ;                              { 评价值=原完工时间-操作期间最好完工时间 }
           end ;

type hls_record = record                                   // 评价高层个体时使用
           arr : array [ 1.. hls_length ] of integer ;     { 一个个算子 }
           op : array [ 1.. max_i ] of integer ;
           eva : double ;
           end ;

type gant_tu = record                                    //甘特图需要记录的数据 ：
           start : arr ;                                 //每个工序块的开始加工时间
           dura : arr ;                                  //每个工序块的加工时间
           macInfo : arr ;                               //每个工序块的工件号
           jobid : arr ;                                 //每个工序块的工序号
           oper : arr ;
           end ;

type double_arr = array [ 1 .. hls_popsize ] of double ; //给高层个体赋值的序列

type SBZM_pop1 = array [ 1 .. initially_generated_individuals ] of SBZM ;             //初始生成x个体

type hls_indiv_pop = array [ 1 .. hls_popsize ] of hls_indiv ;     //高层操作序种群

type hls_indiv_contrast_pop = array [ 1 .. hls_popsize + hls_length ] of hls_indiv ;     //高层操作序种群

type hls_record_pop = array [ 1 .. hls_length ] of hls_record ;    //高层对底层进行操作时，记录高层个体一个一个算子的评价值

type real_three_arr = array [ 1 .. MC_length ] of array [ 1 .. MC_width_and_hight ] of array [ 1 .. MC_width_and_hight ] of real ;  //相似块与概率矩阵



type arr_max = array [ 1 .. 50000 ] of Double ;

/////////////////EDA的变量定义////////////////////////////////////////

type real_two_arr = array [ 1 .. hls_length ] of array [ 1 .. MC_width_and_hight ] of  real ;  //EDA概率矩阵

/////////////////////////////////////////////////////////////////////

  procedure main ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                     source1_match : two_arr4 ;
                                       setup1_time : two_arr1 ;
                          M_mt_arr1 , P_mt_arr1 : three_arr ) ;

  procedure HHEDA ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                     source1_match : two_arr4 ;
                                       setup1_time : two_arr1 ;
                          M_mt_arr1 , P_mt_arr1 : three_arr ) ;

  procedure QHH ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                     source1_match : two_arr4 ;
                                       setup1_time : two_arr1 ;
                          M_mt_arr1 , P_mt_arr1 : three_arr ) ;

  procedure HHMEDA_compare1 ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                     source1_match : two_arr4 ;
                                       setup1_time : two_arr1 ;
                          M_mt_arr1 , P_mt_arr1 : three_arr ) ;




var
  litera_num : Integer ;           //迭代代数
  best_fitness_arr : arr_max ;    //每一代最优加工时间记录

  canshu_arr : array [ 1 .. 30 ] of Double ;

  SBZM_pop , ini_eva_pop ,  pre_pop : SBZM_pop1 ;

  hls_contrast_pop , hls_pop : hls_indiv_contrast_pop ;    // hls_pop ：相似块矩阵对其进行采样
                                                           // hls_contrast_pop ：每生成一代新高层种群，在hls_contrast_pop内进行排序寻优

  best_contrast_ind ,pai0 : SBZM ;     //pai0:初代最优，对其进行一代代改善
                                       //best_contrast_ind ：每代改善完的适配值

  count_xu , start_time , dura_time : arr ;   //甘特图所需数据

  /////////////////机器资源图所需数据///////////////////////////////////////////

  Tester_matrix , handler_matrix : array [ 1 .. 4 ] of array [ 1 .. 240 ] of  integer ;

  Acces_matrix : array [ 1 .. 5 ] of array [ 1 .. 240 ] of  integer ;

  sourcetime : Integer ;

  ////////////////////////////////////////////////////////////////////////////////

  Pro_model_MC : real_three_arr  ;     //概率模型矩阵立方体
  similar_blocks_MC : real_three_arr ;  //相似块矩阵立方体

  point_value : real ;                 //概率指针

  //////////////////EDA变量定义/////////////////////////////////

  EDA_MC : real_two_arr ;

  //////////////////////////////////////////////////////////////





  /////////////////QHH变量定义/////////////////////////////////

  QHH_state , QHH_action , QHH_G , QHH_gcur , QHH_E , QHH_EP : Integer ;

  QHH_pt , QHH_C_init , QHH_C_EP , QHH_greedy , QHH_discount , QHH_learn_rate ,
  QHH_rein_sig ,QHH_argmax_v : Double ;

  pai_b , pai_c : SBZM ;

  QHH_qtable : array [ 1 .. 3 ] of array [ 1 .. QHH_OP ] of double ;

  //////////////////////////////////////////////////////////////////

  time1,time3: cardinal;

implementation



procedure calc_fitness ( temp_o_arr , temp_test1 , temp_hand1 , temp_accessory1 : arr ;
                                                          temp_source1_match : two_arr4 ;
                                                            temp_setup1_time : two_arr1 ;      //机器序列寻优
                                            temp_M_mt_arr1 , temp_P_mt_arr1 : three_arr ;
                                                                   var  o11_arr : arr ;
                                                                   var  m11_arr : arr ;
                                                                   var  pri_arr : arr ;
                                                                  var fitness : Double ) ;

var
  i , temper : Integer ;
  temp_o , trans_o_arr: arr ;
  t_set : set of 1 .. 255 ;
  flag : Boolean ;

  M_arr , P_arr , job_count_arr ,  M_select , gant_count : arr ;
  EC_arr , FC_arr : S_arr ;
  job_M : two_arr2 ;
  insert_flag , op_flag :Boolean ;
  Q , j , k , x , E , Ek , EP , F , Fk , FP , a , n , m , y : Integer ;
  Machine_pro_pos , Machine_pro_pos1 , machine_f_pos , mac_job_pos : Pro_pos ;
  job_Ctime , job_ftime  : three_arr ;
  gant_MT , gant_MT1 : three_arr1 ;
  gant_MT_pos , gant_MT1_pos : two_arr1 ;
  loop : arr_loop ;



  OMP : Indiv1 ;
  t_acces : array [ 1 .. 4 ] of Integer ;
  t_test , t_hand : array [ 1 .. 3 ] of Integer ;
  pr_set : set of 1 .. 36 ;
  m_final : array [ 1 .. 36 ] of Integer ;
  count_op  : array [ 1 .. 100 ] of Integer ;
  count_ma : array [ 1 .. 36 ] of Integer ;
  op_ctime : two_arr2 ;
  op_pos : array [ 1 .. 100 ] of Integer ;
  min_mac , min_mac1 , min_mac2 , min_mac3 : Integer ;
  min_time , min_time1 , min_time2 , min_time3 : Integer ;
  young_arr : two_arr3 ;
  source_count : array [ 1 .. 215 ] of array [ 1 .. 3 ] of array [ 1 .. 4 ] of Integer ;

begin

  fillchar ( job_count_arr , sizeof ( arr ) , 0 ) ;     //记录工序

  fillchar ( job_M , sizeof ( two_arr2 ) , 0 ) ;        //记录工件在哪台机器加工

  fillchar ( job_Ctime , sizeof ( three_arr ) , 0 ) ;   //记录每个工序的开始完工时间

  fillchar ( gant_MT , sizeof ( three_arr1 ) , 0 ) ;

  fillchar ( gant_MT_pos , sizeof ( two_arr1 ) , 0 ) ;        //  用来转化前插解码后甘特图

  fillchar ( gant_MT1_pos , sizeof ( two_arr1 ) , 0 ) ;

  for i := 1 to 36 do

  begin

    fillchar ( machine_pro_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;
                                                                          //machine_pro_pos：记录机器上每个位置的开始结束时间
    fillchar ( machine_pro_pos1 [ i ] , sizeof ( two_arr3 ) , 0 ) ;


  end ;

  for i := 1 to max_i do

  begin

    job_count_arr [ temp_o_arr [ i ] ] := job_count_arr [ temp_o_arr [ i ] ] + 1 ;

    fillchar ( EC_arr , sizeof ( S_arr ) , 0 ) ;                        //EC: 有前插空间时，工序可选择的机器集 ，从中挑选前插空间最小的机器

    fillchar ( FC_arr , sizeof ( S_arr ) , 0 ) ;                        //FC : 无前插机器可选时，挑选加工完成时间最短的机器

    insert_flag := False;

    op_flag := False;

    if job_count_arr [ temp_o_arr [ i ] ] = 1 then                              //如果是工件的第一个工序

    begin

      j := 1 ;

      repeat

        k := 1 ;

        repeat

          if k = 1 then                                                         //轮到找机器第一个工序前有没有空

          begin

            if  temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
                <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]   //如果有前插空间

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k  , 1 ] ;   //②比较插入空间大小

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;              //记录拥有最小前插空间的机器

                EP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;             //记录加工时间

                Ek := k ;                                                       //记录前插的位置

              end ;

            end ;

            if  temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
                > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //没有前插空间

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 2 ] = 0 then    //机器上没有工件（轮到了机器上最后一个工序了）

              begin

                op_flag := True ;

                FC_arr [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;           //记录加工时间

                  Fk := k ;                                                     //记录位置

                end;

              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 2 ] > 0 then     //机器上有工件（还没轮到机器上最后一个工序）

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

          if k > 1 then                                                         //不是机器上加工第一个工序

          begin

            if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]

               + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]

               <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //有前插空间

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k  , 1 ]    //②比较插入空间大小
                           - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] ;

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then                              //选插入空间最小的

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;              //记录拥有最小前插空间的机器

                Ep := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;             //记录加工时间

                EK := k ;                                                       //记录前插的位置

              end ;

            end ;

            if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
               > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //没有前插空间

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]      //机器上没有工件（轮到了机器上最后一个工序了）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] < 0

              then

              begin

                op_flag := True ;

                FC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]
                                + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                  Fk := k ;

                end;
              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //机器上有工件（还没轮到机器上最后一个工序）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] >= 0

              then

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

        until ( insert_flag = True ) or ( op_flag = True ) ;

        j := j + 1 ;

      until ( temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] = 0 ) or ( j = 4 );

      if insert_flag = true then                                                  //有前插空间时选择机器

      begin

        M_arr [ i ] := E ;

        P_arr [ i ] := EP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := E ;

        if Ek = 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;





          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ] := 0 ;                         //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ] := EP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := 0 ;                //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

        if Ek > 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] + Ep ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          :=  machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

      end ;

      if insert_flag = false then                                                 //无前插空间时选择机器

      begin

        M_arr [ i ] := F ;

        P_arr [ i ] := FP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := F ;

        if Fk = 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ] := 0 ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ] := FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := 0 ;                  //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

        if Fk > 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ;                        //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ]  ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

      end ;

    end ;

    if job_count_arr [ temp_o_arr [ i ] ] > 1 then                              //如果不是工件的第一个工序

    begin

      j := 1 ;

      repeat

        k := 1 ;

        repeat

          if K = 1 then                                                         //机器上第一个加工

          begin

            if job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //有前插空间
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

               <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ]
              := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k  , 1 ] ;   //②比较插入空间大小

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;              //记录拥有最小前插空间的机器

                EP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;             //记录加工时间

                Ek := k ;                                                       //记录前插的位置

              end;
            end;

            if  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //没有前插空间
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

               > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 2 ] = 0 then    //到了机器上最后一个工序

              begin

                op_flag := True ;

                FC_arr [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;           //记录加工时间

                  Fk := k ;                                                     //记录位置

                end;

              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 2 ] > 0 then     //机器上有工件（还没轮到机器上最后一个工序）

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

          if K > 1 then                                                         //机器上不是第一个加工

          begin

            if  Max ( job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //有前插空间
                      + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ] ,
                      machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] )
                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

                <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin
              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k  , 1 ]    //②比较插入空间大小
                           - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] ;

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then                              //选插入空间最小的

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;        //记录拥有最小前插空间的机器

                Ep := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;       //记录加工时间

                EK := k ;                                                       //记录前插的位置

              end ;
            end;

            if   Max ( job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //没有前插空间
                      + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ] ,
                      machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] )
                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

                > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]      //机器上没有工件（轮到了机器上最后一个工序了）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] < 0

              then

              begin

                op_flag := True ;

                FC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ]
                                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                  Fk := k ;

                end;
              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]    //机器上有工件（还没轮到机器上最后一个工序）

                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] >= 0

              then

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

        until  (insert_flag = True) or (op_flag = True) ;

        j := j + 1 ;

      until ( temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] = 0 ) or  ( j = 4 ) ;

      if insert_flag = True then                                                //可选机器内有可前插的
      begin

        M_arr [ i ] := E ;

        P_arr [ i ] := EP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := E ;

        if Ek = 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;            //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ] :=  machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + EP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

        if Ek > 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          :=  Max( machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          :=  machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

      end ;

      if insert_flag = false then                                                //可选机器内有可前插的

      begin

        M_arr [ i ] := F ;

        P_arr [ i ] := FP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := F ;

        if Fk = 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ] := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                  //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

        if Fk > 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] )  ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //⑤更新工件i第n道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end ;

      end ;

    end ;

  end ;

  ////////////////////////////////////////////////////////////////////////

  n := 1 ;

  m := 1 ;

  fillchar ( gant_count , sizeof ( arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    loop [ i ] := 1 ;

  end;

  repeat

    if gant_MT [ m , loop [ m ] , 1 ] = 0 then

    begin

      m := m + 1 ;

    end

    else

    begin

      if gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] = gant_MT [ m , loop [ m ] , 2 ] - 1 then

      begin

        gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] := gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] + 1 ;

        OMP.op_arr [ n ] := gant_MT [ m , loop [ m ] , 1 ] ;

        OMP.ma_arr [ n ] := m ;

        OMP.pr_arr [ n ] := temp_p_mt_arr1 [ gant_MT [ m , loop [ m ] , 1 ] ,
                                             gant_MT [ m , loop [ m ] , 2 ] ,
                                             gant_MT_pos [ m , loop [ m ] ] ] ;

        loop [ m ] := loop [ m ] + 1 ;

        m := m + 1 ;

        n := n + 1 ;

      end

      else

      begin

        m := m + 1 ;

      end;

    end;

    if m >= 37 then

    begin

      m := 1 ;

    end;

  until n > max_i ;

  for i := 1 to 3 do

  begin

    t_test [ i ] := temp_test1 [ i ] ;

    t_hand [ i ] := temp_hand1 [ i ] ;

  end ;

  for i := 1 to 4 do

  begin

    t_acces [ i ] := temp_accessory1 [ i ] ;

  end ;

  pr_set := [ ] ;

  fillchar ( count_op , sizeof ( count_op ) , 0 ) ;

  fillchar ( count_ma , sizeof ( count_ma ) , 0 ) ;

  fillchar ( op_ctime , sizeof ( op_ctime ) , 0 ) ;

  fillchar ( m_final , sizeof ( m_final ) , 0 ) ;

  fillchar ( op_pos , sizeof ( op_pos ) , 0 ) ;

  fillchar ( job_ftime , sizeof ( three_arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    fillchar ( machine_f_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

    fillchar ( mac_job_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

  end ;

  for i := 1 to 215 do

  begin

    for j := 1 to 3 do

    begin

      for k := 1 to 4 do

      begin
        source_count [ i , j , k ] := 0 ;

      end;

    end;
    
  end;

////////////////////////////////////////////////                                加上资源约束计算完工时间

  for i := 1 to max_i do

  begin

    count_op [ OMP.op_arr [ i ] ] := count_op [ OMP.op_arr [ i ] ] + 1 ;        //工序+1

    count_ma [ OMP.ma_arr [ i ] ] := count_ma [ OMP.ma_arr [ i ] ] + 1 ;        //机器加工位置+1

    min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

    fillchar ( young_arr , sizeof ( two_arr3 ) , 0 ) ;

    if count_ma [ OMP.ma_arr [ i ] ] = 1 then                                   //机器第一个位置加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )     //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] ;               //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] := 0 ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := OMP.pr_arr [ i ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + min_time ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                          //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := OMP.op_arr [ i ] ;

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
           := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;


          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器第一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]                 //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                                  //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;             //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;             //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

           y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

           t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
           := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

           t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
           := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

           t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
           := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

           if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

           begin

             pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

           end;

        end ;

      end ;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    end
    else                                   //不是在机器第一位工序加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + m_final [ OMP.ma_arr [ i ] ] ;          //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ]
           + max ( min_time , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end
      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //不是机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ,  m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end

      else

      if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器非一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] :=
          OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin


                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end ;

      if count_op [ OMP.op_arr [ i ] ] = 4 then                            //机器非一个位置加工，工件第四个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] :=
          OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin


                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end ;

    end ;

    source_count [ i , 1 , 1 ] := t_test [ 1 ] ;

    source_count [ i , 1 , 2 ] := t_test [ 2 ] ;

    source_count [ i , 1 , 3 ] := t_test [ 3 ] ;

    source_count [ i , 2 , 1 ] := t_hand [ 1 ] ;

    source_count [ i , 2 , 2 ] := t_hand [ 2 ] ;

    source_count [ i , 2 , 3 ] := t_hand [ 3 ] ;

    source_count [ i , 3 , 1 ] := t_acces [ 1 ] ;

    source_count [ i , 3 , 2 ] := t_acces [ 2 ] ;

    source_count [ i , 3 , 3 ] := t_acces [ 3 ] ;

    source_count [ i , 3 , 4 ] := t_acces [ 4 ] ;

  end ;



//////////////////////////////////////////////////////////////////////////////////////

  fitness := m_final [ 1 ] ;

  for i := 2 to 36 do

  begin

    if fitness < m_final [ i ] then

    begin

      fitness := m_final [ i ] ;

    end;

  end;

  for i := 1 to max_i do

  begin

    o11_arr [ i ] := OMP.op_arr [ i ] ;

    m11_arr [ i ] := OMP.ma_arr [ i ] ;

    pri_arr [ i ] := temp_o_arr [ i ] ;

  end;

end ;

procedure calc_fit_and_data_stored ( temp_o_arr , temp_test1 , temp_hand1 , temp_accessory1 : arr ;
                                                          temp_source1_match : two_arr4 ;
                                                            temp_setup1_time : two_arr1 ;
                                            temp_M_mt_arr1 , temp_P_mt_arr1 : three_arr ;
                                                                   var  o11_arr : arr ;
                                                                   var  m11_arr : arr ;
                                                                   var  pri_arr : arr ;
                                                                  var fitness : Double ) ;


var
  i , temper : Integer ;
  temp_o , trans_o_arr: arr ;
  t_set : set of 1 .. 255 ;
  flag : Boolean ;

  M_arr , P_arr , job_count_arr ,  M_select , gant_count : arr ;
  EC_arr , FC_arr : S_arr ;
  job_M : two_arr2 ;
  insert_flag , op_flag :Boolean ;
  Q , j , k , x , E , Ek , EP , F , Fk , FP , a , n , m , y : Integer ;
  Machine_pro_pos , Machine_pro_pos1 , machine_f_pos , mac_job_pos : Pro_pos ;
  job_Ctime , job_ftime  : three_arr ;
  gant_MT , gant_MT1 : three_arr1 ;
  gant_MT_pos , gant_MT1_pos : two_arr1 ;
  loop : arr_loop ;

  OMP : Indiv1 ;
  t_acces : array [ 1 .. 4 ] of Integer ;
  t_test , t_hand : array [ 1 .. 3 ] of Integer ;
  pr_set : set of 1 .. 36 ;
  m_final : array [ 1 .. 36 ] of Integer ;
  count_op  : array [ 1 .. 100 ] of Integer ;
  count_ma : array [ 1 .. 36 ] of Integer ;
  op_ctime : two_arr2 ;
  op_pos : array [ 1 .. 100 ] of Integer ;
  min_mac , min_mac1 , min_mac2 , min_mac3 : Integer ;
  min_time , min_time1 , min_time2 , min_time3 : Integer ;
  young_arr : two_arr3 ;
  source_count : array [ 1 .. 215 ] of array [ 1 .. 3 ] of array [ 1 .. 4 ] of Integer ;

begin

  fillchar ( job_count_arr , sizeof ( arr ) , 0 ) ;

  fillchar ( job_M , sizeof ( two_arr2 ) , 0 ) ;

  fillchar ( job_Ctime , sizeof ( three_arr ) , 0 ) ;

  fillchar ( gant_MT , sizeof ( three_arr1 ) , 0 ) ;

  fillchar ( gant_MT_pos , sizeof ( two_arr1 ) , 0 ) ;

  fillchar ( gant_MT1_pos , sizeof ( two_arr1 ) , 0 ) ;

  /////////////////////////////////////////////////////////////////

  fillchar ( count_xu , sizeof ( arr ) , 0 ) ;

  fillchar ( start_time , sizeof ( arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    fillchar ( machine_pro_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

    fillchar ( machine_pro_pos1 [ i ] , sizeof ( two_arr3 ) , 0 ) ;


  end ;

  for i := 1 to max_i do

  begin

    job_count_arr [ temp_o_arr [ i ] ] := job_count_arr [ temp_o_arr [ i ] ] + 1 ;

    fillchar ( EC_arr , sizeof ( S_arr ) , 0 ) ;

    fillchar ( FC_arr , sizeof ( S_arr ) , 0 ) ;

    insert_flag := False;

    op_flag := False;

    if job_count_arr [ temp_o_arr [ i ] ] = 1 then                              //如果是工件的第一个工序

    begin

      j := 1 ;

      repeat

        k := 1 ;

        repeat

          if k = 1 then                                                         //轮到找机器第一个工序前有没有空

          begin

            if  temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
                <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]   //如果有前插空间

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k  , 1 ] ;   //②比较插入空间大小

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;              //记录拥有最小前插空间的机器

                EP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;             //记录加工时间

                Ek := k ;                                                       //记录前插的位置

              end ;

            end ;

            if  temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
                > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //没有前插空间

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 2 ] = 0 then    //机器上没有工件（轮到了机器上最后一个工序了）

              begin

                op_flag := True ;

                FC_arr [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;           //记录加工时间

                  Fk := k ;                                                     //记录位置

                end;

              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 2 ] > 0 then     //机器上有工件（还没轮到机器上最后一个工序）

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

          if k > 1 then                                                         //不是机器上加工第一个工序

          begin

            if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]

               + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]

               <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //有前插空间

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k  , 1 ]    //②比较插入空间大小
                           - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] ;

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then                              //选插入空间最小的

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;              //记录拥有最小前插空间的机器

                Ep := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;             //记录加工时间

                EK := k ;                                                       //记录前插的位置

              end ;

            end ;

            if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ]
               > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //没有前插空间

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]      //机器上没有工件（轮到了机器上最后一个工序了）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] < 0

              then

              begin

                op_flag := True ;

                FC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ]
                                + temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] ;

                  Fk := k ;

                end;
              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k , 1 ]    //机器上有工件（还没轮到机器上最后一个工序）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] , k - 1 , 2 ] >= 0

              then

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

        until ( insert_flag = True ) or ( op_flag = True ) ;

        j := j + 1 ;

      until ( temp_M_mt_arr1 [ temp_o_arr [ i ] , 1 , j ] = 0 ) or ( j = 4 );

      if insert_flag = true then                                                  //有前插空间时选择机器

      begin

        M_arr [ i ] := E ;

        P_arr [ i ] := EP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := E ;

        if Ek = 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;





          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ] := 0 ;                         //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ] := EP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := 0 ;                //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

        if Ek > 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] + Ep ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          :=  machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

      end ;

      if insert_flag = false then                                                 //无前插空间时选择机器

      begin

        M_arr [ i ] := F ;

        P_arr [ i ] := FP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := F ;

        if Fk = 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ] := 0 ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ] := FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := 0 ;                  //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

        if Fk > 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ;                        //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ]  ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

      end ;

    end ;

    if job_count_arr [ temp_o_arr [ i ] ] > 1 then                              //如果不是工件的第一个工序

    begin

      j := 1 ;

      repeat

        k := 1 ;

        repeat

          if K = 1 then                                                         //机器上第一个加工

          begin

            if job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //有前插空间
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

               <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ]
              := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k  , 1 ] ;   //②比较插入空间大小

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;              //记录拥有最小前插空间的机器

                EP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;             //记录加工时间

                Ek := k ;                                                       //记录前插的位置

              end;
            end;

            if  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //没有前插空间
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ]
               + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

               > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 2 ] = 0 then    //到了机器上最后一个工序

              begin

                op_flag := True ;

                FC_arr [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;           //记录加工时间

                  Fk := k ;                                                     //记录位置

                end;

              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 2 ] > 0 then     //机器上有工件（还没轮到机器上最后一个工序）

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

          if K > 1 then                                                         //机器上不是第一个加工

          begin

            if  Max ( job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //有前插空间
                      + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ] ,
                      machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] )
                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

                <= machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin
              insert_flag := True ;                                             //①改变前插标识

              EC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k  , 1 ]    //②比较插入空间大小
                           - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] ;

              if EC_arr [ 1 ] = 0 then

              begin

                EC_arr [ 1 ] := EC_arr [ j ] ;

              end;

              if EC_arr [ j ] <= EC_arr [ 1 ]  then                              //选插入空间最小的

              begin

                Q := j ;

                EC_arr [ 1 ] := EC_arr [ j ] ;

                E := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;        //记录拥有最小前插空间的机器

                Ep := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;       //记录加工时间

                EK := k ;                                                       //记录前插的位置

              end ;
            end;

            if   Max ( job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]        //没有前插空间
                      + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ] ,
                      machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] )
                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ]

                > machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]

            then

            begin

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]      //机器上没有工件（轮到了机器上最后一个工序了）
                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] < 0

              then

              begin

                op_flag := True ;

                FC_arr [ j ] := machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ]
                                + temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                if FC_arr [ j ] <= FC_arr [ 1 ] then

                begin

                  Q := j ;

                  FC_arr [ 1 ] := FC_arr [ j ] ;

                  F := temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;            //选定拥有最小完工时间的机器

                  FP := temp_P_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] ;

                  Fk := k ;

                end;
              end ;

              if machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k , 1 ]    //机器上有工件（还没轮到机器上最后一个工序）

                 - machine_pro_pos [ temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] , k - 1 , 2 ] >= 0

              then

              begin

                k := k + 1 ;

              end ;

            end ;

          end ;

        until  (insert_flag = True) or (op_flag = True) ;

        j := j + 1 ;

      until ( temp_M_mt_arr1 [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , j ] = 0 ) or  ( j = 4 ) ;

      if insert_flag = True then                                                //可选机器内有可前插的
      begin

        M_arr [ i ] := E ;

        P_arr [ i ] := EP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := E ;

        if Ek = 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;            //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ] :=  machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + EP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

        if Ek > 1 then

        begin

          for a := 1 to 36 do

          begin

            machine_pro_pos1 [ M_arr [ i ] , a , 1 ]
            := machine_pro_pos [ M_arr [ i ] , a , 1 ] ;

            machine_pro_pos1 [ M_arr [ i ] , a , 2 ]
            := machine_pro_pos [ M_arr [ i ] , a , 2 ] ;

            gant_MT1 [ M_arr [ i ] , a , 1 ]
            := gant_MT [ M_arr [ i ] , a , 1 ] ;

            gant_MT1 [ M_arr [ i ] , a , 2 ]
            := gant_MT [ M_arr [ i ] , a , 2 ] ;

            gant_MT1_pos [ M_arr [ i ] , a ]
            := gant_MT_pos [ M_arr [ i ] , a ] ;

          end;

          for x := Ek to 35 do

          begin                                                                   //③确定机器后改变机器从第k个开始顺延的加工工件的开始结束时间

            machine_pro_pos [ M_arr [ i ] , x + 1 , 1 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 1 ] ;

            machine_pro_pos [ M_arr [ i ] , x + 1 , 2 ]
            := machine_pro_pos1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 1 ]
            := gant_MT1 [ M_arr [ i ] , x , 1 ] ;

            gant_MT [ M_arr [ i ] , x + 1 , 2 ]
            := gant_MT1 [ M_arr [ i ] , x , 2 ] ;

            gant_MT_pos [ M_arr [ i ] , x + 1 ]
            := gant_MT1_pos [ M_arr [ i ] , x ] ;

          end;

          machine_pro_pos [ M_arr [ i ] , Ek , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Ek , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          :=  Max( machine_pro_pos [ M_arr [ i ] , Ek - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          :=  machine_pro_pos [ M_arr [ i ] , Ek , 1 ] + Ep ;

          gant_MT [ M_arr [ i ] , Ek , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Ek , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Ek ] := Q ;

        end;

      end ;

      if insert_flag = false then                                                //可选机器内有可前插的

      begin

        M_arr [ i ] := F ;

        P_arr [ i ] := FP ;

        job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] ] := F ;

        if Fk = 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ] := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ] := job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
               + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ;                  //⑤更新工件i第一道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ] := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end;

        if Fk > 1 then

        begin

          machine_pro_pos [ M_arr [ i ] , Fk , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] )  ;                       //④更新机器j第k个的加工工序开始结束时间

          machine_pro_pos [ M_arr [ i ] , Fk , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 1 ]
          := Max( machine_pro_pos [ M_arr [ i ] , Fk - 1 , 2 ] ,
                  job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 , 2 ]
                  + temp_setup1_time [ job_M [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] - 1 ] ,
                                    M_arr [ i ] ] ) ;                       //⑤更新工件i第n道工序的加工完成时间

          job_Ctime [ temp_o_arr [ i ] , job_count_arr [ temp_o_arr [ i ] ] , 2 ]
          := machine_pro_pos [ M_arr [ i ] , Fk , 1 ] + FP ;

          gant_MT [ M_arr [ i ] , Fk , 1 ] := temp_o_arr [ i ] ;

          gant_MT [ M_arr [ i ] , Fk , 2 ] := job_count_arr [ temp_o_arr [ i ] ] ;

          gant_MT_pos [ M_arr [ i ] , Fk ] := Q ;

        end ;

      end ;

    end ;

  end ;

  ////////////////////////////////////////////////////////////////////////

  n := 1 ;

  m := 1 ;

  fillchar ( gant_count , sizeof ( arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    loop [ i ] := 1 ;

  end;

  repeat

    if gant_MT [ m , loop [ m ] , 1 ] = 0 then

    begin

      m := m + 1 ;

    end

    else

    begin

      if gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] = gant_MT [ m , loop [ m ] , 2 ] - 1 then

      begin

        gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] := gant_count [ gant_MT [ m , loop [ m ] , 1 ] ] + 1 ;

        OMP.op_arr [ n ] := gant_MT [ m , loop [ m ] , 1 ] ;

        OMP.ma_arr [ n ] := m ;

        OMP.pr_arr [ n ] := temp_p_mt_arr1 [ gant_MT [ m , loop [ m ] , 1 ] ,
                                             gant_MT [ m , loop [ m ] , 2 ] ,
                                             gant_MT_pos [ m , loop [ m ] ] ] ;

        loop [ m ] := loop [ m ] + 1 ;

        m := m + 1 ;

        n := n + 1 ;

      end

      else

      begin

        m := m + 1 ;

      end;

    end;

    if m >= 37 then

    begin

      m := 1 ;

    end;

  until n > max_i ;

  for i := 1 to 3 do

  begin

    t_test [ i ] := temp_test1 [ i ] ;

    t_hand [ i ] := temp_hand1 [ i ] ;

  end ;

  for i := 1 to 4 do

  begin

    t_acces [ i ] := temp_accessory1 [ i ] ;

  end ;

  pr_set := [ ] ;

  fillchar ( count_op , sizeof ( count_op ) , 0 ) ;

  fillchar ( count_ma , sizeof ( count_ma ) , 0 ) ;

  fillchar ( op_ctime , sizeof ( op_ctime ) , 0 ) ;

  fillchar ( m_final , sizeof ( m_final ) , 0 ) ;

  fillchar ( op_pos , sizeof ( op_pos ) , 0 ) ;

  fillchar ( job_ftime , sizeof ( three_arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    fillchar ( machine_f_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

    fillchar ( mac_job_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

  end ;

  for i := 1 to 215 do

  begin

    for j := 1 to 3 do

    begin

      for k := 1 to 4 do

      begin
        source_count [ i , j , k ] := 0 ;

      end;

    end;
    
  end;

 /////////////////////初始化机器资源数组/////////////////////////////////////

  for i := 1 to 240 do

  begin

    Tester_matrix [ 1 , i ] := i - 1 ;

    Tester_matrix [ 2 , i ] := -1 ;

    Tester_matrix [ 3 , i ] := -1 ;

    Tester_matrix [ 4 , i ] := -1 ;

    handler_matrix [ 1 , i ] := i - 1 ;

    handler_matrix [ 2 , i ] := -1 ;

    handler_matrix [ 3 , i ] := -1 ;

    handler_matrix [ 4 , i ] := -1 ;

    Acces_matrix [ 1 , i ] := i - 1 ;

    Acces_matrix [ 2 , i ] := -1 ;

    Acces_matrix [ 3 , i ] := -1 ;

    Acces_matrix [ 4 , i ] := -1 ;

    Acces_matrix [ 5 , i ] := -1 ;

  end;

////////////////////////////////////////////////                                加上资源约束计算完工时间

  for i := 1 to max_i do

  begin

    count_op [ OMP.op_arr [ i ] ] := count_op [ OMP.op_arr [ i ] ] + 1 ;        //工序+1

    count_ma [ OMP.ma_arr [ i ] ] := count_ma [ OMP.ma_arr [ i ] ] + 1 ;        //机器加工位置+1

    min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

    fillchar ( young_arr , sizeof ( two_arr3 ) , 0 ) ;

    if count_ma [ OMP.ma_arr [ i ] ] = 1 then                                   //机器第一个位置加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )     //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] ;               //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] := 0 ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := OMP.pr_arr [ i ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + min_time ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                          //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := OMP.op_arr [ i ] ;

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
           := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;


          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器第一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]                 //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                                  //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;             //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;             //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

           y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

           t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
           := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

           t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
           := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

           t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
           := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

           if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

           begin

             pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

           end;

        end ;

      end ;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    end
    else                                   //不是在机器第一位工序加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + m_final [ OMP.ma_arr [ i ] ] ;          //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ]
           + max ( min_time , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end
      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //不是机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ,  m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end

      else

      if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器非一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] :=
          OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin


                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;
                     
                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end ;

    end ;

    source_count [ i , 1 , 1 ] := t_test [ 1 ] ;

    source_count [ i , 1 , 2 ] := t_test [ 2 ] ;

    source_count [ i , 1 , 3 ] := t_test [ 3 ] ;

    source_count [ i , 2 , 1 ] := t_hand [ 1 ] ;

    source_count [ i , 2 , 2 ] := t_hand [ 2 ] ;

    source_count [ i , 2 , 3 ] := t_hand [ 3 ] ;

    source_count [ i , 3 , 1 ] := t_acces [ 1 ] ;

    source_count [ i , 3 , 2 ] := t_acces [ 2 ] ;

    source_count [ i , 3 , 3 ] := t_acces [ 3 ] ;

    source_count [ i , 3 , 4 ] := t_acces [ 4 ] ;

    start_time [ i ] := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;

    count_xu [ i ] := count_op [ OMP.op_arr [ i ] ] ;

  //////////////////机器资源数据///////////////////////////////////////

    tester_matrix [ 2 , start_time [ i ] + 1 ] := t_test [ 1 ] ;

    tester_matrix [ 3 , start_time [ i ] + 1 ] := t_test [ 2 ] ;

    tester_matrix [ 4 , start_time [ i ] + 1 ] := t_test [ 3 ] ;

    handler_matrix [ 2 , start_time [ i ] + 1 ] := t_hand [ 1 ] ;

    handler_matrix [ 3 , start_time [ i ] + 1 ] := t_hand [ 2 ] ;

    handler_matrix [ 4 , start_time [ i ] + 1 ] := t_hand [ 3 ] ;

    acces_matrix [ 2 , start_time [ i ] + 1 ] := t_acces [ 1 ] ;

    acces_matrix [ 3 , start_time [ i ] + 1 ] := t_acces [ 2 ] ;

    acces_matrix [ 4 , start_time [ i ] + 1 ] := t_acces [ 3 ] ;

    acces_matrix [ 5 , start_time [ i ] + 1 ] := t_acces [ 4 ] ;

  end ;

 ///////////////////////////////////////////////////////////////////////

//////////////甘特图,机器资源图数据/////////////////////////////////////////////////

  fitness := m_final [ 1 ] ;

  sourcetime := m_final [ 1 ] ;

  for i := 2 to 36 do

  begin

    if fitness < m_final [ i ] then

    begin

      fitness := m_final [ i ] ;

      sourcetime := m_final [ i ] ;

    end;

  end;

  for i := 1 to max_i do

  begin

    o11_arr [ i ] := OMP.op_arr [ i ] ;

    m11_arr [ i ] := OMP.ma_arr [ i ] ;

    pri_arr [ i ] := temp_o_arr [ i ] ;

    dura_time [ i ] := OMP.pr_arr [ i ] ;

  end;

//////////////////机器资源数据///////////////////////////////////////

  sourcetime := sourcetime + 1 ;
  
  for i := 1 to sourcetime do

  begin

    if ( tester_matrix [ 2 , i ] = -1 ) then

    begin

      tester_matrix [ 2 , i ] := tester_matrix [ 2 , i - 1 ] ;

    end;

    if ( tester_matrix [ 3 , i ] = -1 ) then

    begin

      tester_matrix [ 3 , i ] := tester_matrix [ 3 , i - 1 ] ;

    end;

    if ( tester_matrix [ 4 , i ] = -1 ) then

    begin

      tester_matrix [ 4 , i ] := tester_matrix [ 4 , i - 1 ] ;

    end;

    if ( handler_matrix [ 2 , i ] = -1 ) then

    begin

      handler_matrix [ 2 , i ] := handler_matrix [ 2 , i - 1 ] ;

    end;

    if ( handler_matrix [ 3 , i ] = -1 ) then

    begin

      handler_matrix [ 3 , i ] := handler_matrix [ 3 , i - 1 ] ;

    end;

    if ( handler_matrix [ 4 , i ] = -1 ) then

    begin

      handler_matrix [ 4 , i ] := handler_matrix [ 4 , i - 1 ] ;

    end;

    if ( acces_matrix [ 2 , i ] = -1 ) then

    begin

      acces_matrix [ 2 , i ] := acces_matrix [ 2 , i - 1 ] ;

    end;

    if ( acces_matrix [ 3 , i ] = -1 ) then

    begin

      acces_matrix [ 3 , i ] := acces_matrix [ 3 , i - 1 ] ;

    end;

    if ( acces_matrix [ 4 , i ] = -1 ) then

    begin

      acces_matrix [ 4 , i ] := acces_matrix [ 4 , i - 1 ] ;

    end;

    if ( acces_matrix [ 5 , i ] = -1 ) then

    begin

      acces_matrix [ 5 , i ] := acces_matrix [ 5 , i - 1 ] ;

    end;

  end;

end ;

procedure calc_fitness_compare ( temp_o_arr , temp_test1 , temp_hand1 , temp_accessory1 : arr ;
                                                          temp_source1_match : two_arr4 ;
                                                            temp_setup1_time : two_arr1 ;      //机器序列寻优
                                            temp_M_mt_arr1 , temp_P_mt_arr1 : three_arr ;
                                                                   var  o11_arr : arr ;
                                                                   var  m11_arr : arr ;
                                                                   var  pri_arr : arr ;
                                                                  var fitness : Double ) ;

var
  i , temper : Integer ;
  temp_o , trans_o_arr: arr ;
  t_set : set of 1 .. 255 ;
  flag : Boolean ;

  M_arr , P_arr , job_count_arr ,  M_select , gant_count : arr ;
  EC_arr , FC_arr : S_arr ;
  job_M : two_arr2 ;
  insert_flag , op_flag :Boolean ;
  Q , j , k , x , E , Ek , EP , F , Fk , FP , a , n , m , y : Integer ;
  Machine_pro_pos , Machine_pro_pos1 , machine_f_pos , mac_job_pos : Pro_pos ;
  job_Ctime , job_ftime  : three_arr ;
  gant_MT , gant_MT1 : three_arr1 ;
  gant_MT_pos , gant_MT1_pos : two_arr1 ;
  loop : arr_loop ;

  ///////////compare///////////////

  op_num , job_makespan , mac_makespan , mac , complete , protime: arr ;     //工序，工件完工时间，机器最后一个工件完工时间 ,机器，完工时间，加工时间，
  op_select_mac : two_arr2 ;                       //工序选择加工的机器
  mac_select : Integer ;           //机器选择

  ///////////compare///////////////

  OMP : Indiv1 ;
  t_acces : array [ 1 .. 4 ] of Integer ;
  t_test , t_hand : array [ 1 .. 3 ] of Integer ;
  pr_set : set of 1 .. 36 ;
  m_final : array [ 1 .. 36 ] of Integer ;
  count_op  : array [ 1 .. 100 ] of Integer ;
  count_ma : array [ 1 .. 36 ] of Integer ;
  op_ctime : two_arr2 ;
  op_pos : array [ 1 .. 100 ] of Integer ;
  min_mac , min_mac1 , min_mac2 , min_mac3 : Integer ;
  min_time , min_time1 , min_time2 , min_time3 : Integer ;
  young_arr : two_arr3 ;
  source_count : array [ 1 .. 215 ] of array [ 1 .. 3 ] of array [ 1 .. 4 ] of Integer ;

begin

  fillchar ( op_num , sizeof ( arr ) , 0 ) ;

  fillchar ( job_makespan , sizeof ( arr ) , 0 ) ;

  fillchar ( mac_makespan , sizeof ( arr ) , 0 ) ;

  fillchar ( op_select_mac , sizeof ( two_arr2 ) , 0 ) ;

  for i := 1 to max_i do

  begin

    op_num [ temp_o_arr [ i ] ] := op_num [ temp_o_arr [ i ] ] + 1 ;   //工序+1

    j := 1 ;                                                           //遍历工序所有可选择加工机器

    fillchar ( complete , sizeof ( arr ) , 0 ) ;

    fillchar ( mac , sizeof ( arr ) , 0 ) ;

    repeat

      if ( op_num [ temp_o_arr [ i ] ] = 1 ) then

      begin

        complete [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ]
                                           + mac_makespan [ temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ] ;

        mac [ j ] := temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ;

        protime [ j ] :=  temp_P_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ;

        if ( complete [ j ] < complete [ 1 ] ) then

        begin

          complete [ 1 ] := complete [ j ] ;

          mac [ 1 ] := mac [ j ] ;

          protime [ 1 ] := protime [ j ] ;

        end;

      end;

      if ( op_num [ temp_o_arr [ i ] ] > 1 ) then

      begin

        complete [ j ] := temp_P_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ]
                        + max ( mac_makespan [ temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ] ,
                                job_makespan [ temp_o_arr [ i ] ]
                                + temp_setup1_time [ op_select_mac [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] - 1 ] ,
                                                     temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ] ) ;

        mac [ j ] := temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ;

        protime [ j ] :=  temp_P_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] ;

        if ( complete [ j ] < complete [ 1 ] ) then

        begin

          complete [ 1 ] := complete [ j ] ;

          mac [ 1 ] := mac [ j ] ;

          protime [ 1 ] := protime [ j ] ;

        end;

      end;

      j := j + 1 ;

    until ( temp_M_mt_arr1 [ temp_o_arr [ i ] , op_num [ temp_o_arr [ i ] ] , j ] = 0 ) or ( j = 4 ) ;

    job_makespan [ temp_o_arr [ i ] ] := complete [ 1 ] ;

    mac_makespan [ mac [ 1 ] ] := complete [ 1 ] ;

    OMP.op_arr [ i ] := temp_o_arr [ i ] ;

    OMP.ma_arr [ i ] := mac [ 1 ] ;

    OMP.pr_arr [ i ] := protime [ 1 ] ;

  end ;

  for i := 1 to 3 do

  begin

    t_test [ i ] := temp_test1 [ i ] ;

    t_hand [ i ] := temp_hand1 [ i ] ;

  end ;

  for i := 1 to 4 do

  begin

    t_acces [ i ] := temp_accessory1 [ i ] ;

  end ;

  pr_set := [ ] ;

  fillchar ( count_op , sizeof ( count_op ) , 0 ) ;

  fillchar ( count_ma , sizeof ( count_ma ) , 0 ) ;

  fillchar ( op_ctime , sizeof ( op_ctime ) , 0 ) ;

  fillchar ( m_final , sizeof ( m_final ) , 0 ) ;

  fillchar ( op_pos , sizeof ( op_pos ) , 0 ) ;

  fillchar ( job_ftime , sizeof ( three_arr ) , 0 ) ;

  for i := 1 to 36 do

  begin

    fillchar ( machine_f_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

    fillchar ( mac_job_pos [ i ] , sizeof ( two_arr3 ) , 0 ) ;

  end ;

  for i := 1 to 215 do

  begin

    for j := 1 to 3 do

    begin

      for k := 1 to 4 do

      begin
        source_count [ i , j , k ] := 0 ;

      end;

    end;
    
  end;

////////////////////////////////////////////////                                加上资源约束计算完工时间

  for i := 1 to max_i do

  begin

    count_op [ OMP.op_arr [ i ] ] := count_op [ OMP.op_arr [ i ] ] + 1 ;        //工序+1

    count_ma [ OMP.ma_arr [ i ] ] := count_ma [ OMP.ma_arr [ i ] ] + 1 ;        //机器加工位置+1

    min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

    fillchar ( young_arr , sizeof ( two_arr3 ) , 0 ) ;

    if count_ma [ OMP.ma_arr [ i ] ] = 1 then                                   //机器第一个位置加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )     //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] ;               //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] := 0 ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := OMP.pr_arr [ i ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + min_time ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                          //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := OMP.op_arr [ i ] ;

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
           := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;


          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end;

      end

      else if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器第一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ;                        //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]                 //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                                  //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;             //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;             //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

           job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
           := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

           mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

           y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

           t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
           := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

           t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
           := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

           t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
           := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

           if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

           begin

             pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

           end;

        end ;

      end ;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    end
    else                                   //不是在机器第一位工序加工

    begin

      if count_op [ OMP.op_arr [ i ] ] = 1 then                                 //工件第一个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + m_final [ OMP.ma_arr [ i ] ] ;          //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin

                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ]
           + max ( min_time , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end
      else if count_op [ OMP.op_arr [ i ] ] = 2 then                            //不是机器第一个位置加工，工件第二个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ]
          := OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin



                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) ,  m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end

      else

      if count_op [ OMP.op_arr [ i ] ] = 3 then                            //机器非一个位置加工，工件第三个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] :=
          OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin


                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end ;

      if count_op [ OMP.op_arr [ i ] ] = 4 then                            //机器非一个位置加工，工件第四个工序

      begin

        if   ( t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] > 0 )
         and ( t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] > 0 )    //资源充足时
         and ( t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] > 0 )
        then

        begin

          op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] :=
          OMP.pr_arr [ i ] + max ( op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ]
             + temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ,
               m_final [ OMP.ma_arr [ i ] ] ) ;                                 //完工时间

          m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

          op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] -  OMP.pr_arr [ i ] ;  //

          machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
          := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end

        else

        begin

           y := 1 ;

           min_mac := 0 ; min_mac1 := 0 ; min_mac2 := 0 ; min_mac3 := 0 ;

           min_time := 0 ; min_time1 := 0 ; min_time2 := 0 ; min_time3 := 0 ;

           if t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] = 0

           then

           begin

             min_time1 := 1000 ;

             for j := 1 to 36 do

             begin

               if ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 1 ] = temp_source1_match [ j , 1 ] )

                 then

                 begin



                   if m_final [ j ] < min_time1 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time1 :=  m_final [ j ] ;

                     min_mac1 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           if t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] = 0 then

           begin

             min_time2 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin



                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 2 ] = temp_source1_match [ j , 2 ] )

                 then

                 begin



                   if m_final [ j ] < min_time2 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time2 :=  m_final [ j ] ;

                     min_mac2 := j ;

                   end;

                 end;

               end;

             end;

           end;

           if t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] = 0 then

           begin

             min_time3 := 1000 ;

             for j := 1 to 36 do

             begin

               if  ( j in pr_set )

               then

               begin


                 if ( temp_source1_match [ OMP.ma_arr [ i ] , 3 ] = temp_source1_match [ j , 3 ] )

                 then

                 begin


                   if m_final [ j ] < min_time3 then

                   begin

                     young_arr [ y , 1 ] := j ;

                     young_arr [ y , 2 ] := m_final [ j ] ;

                     y := y + 1 ;

                     min_time3 :=  m_final [ j ] ;

                     min_mac3 := j ;

                   end;

                 end ;

               end;

             end;

           end;

           min_time := min_time1 ;

           min_mac := min_mac1 ;

           if min_time < min_time2 then

           begin

             min_time := min_time2 ;

             min_mac := min_mac2 ;

           end;

           if min_time < min_time3 then

           begin

             min_time := min_time3 ;

             min_mac := min_mac3 ;

           end;

           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] := OMP.pr_arr [ i ] + max ( max ( min_time ,
           op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] - 1 ] +
           temp_setup1_time [ op_pos [ OMP.op_arr [ i ] ] , OMP.ma_arr [ i ] ] ) , m_final [ OMP.ma_arr [ i ] ] ) ;

           m_final [ OMP.ma_arr [ i ] ] := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;   //机器最后一个工件完工时间

           op_pos [ OMP.op_arr [ i ] ] := OMP.ma_arr [ i ] ;                     //此工件加工的机器号

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] - OMP.pr_arr [ i ] ;  //

           machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]         //
           := op_ctime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] ] ;
                                                                                         //更新机器矩阵和工件矩阵
          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 1 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ] ;    //

          job_ftime [ OMP.op_arr [ i ] , count_op [ OMP.op_arr [ i ] ] , 2 ]
          := machine_f_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ] ;    //

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 1 ]
          := OMP.op_arr [ i ] ;

          mac_job_pos [ OMP.ma_arr [ i ] , count_ma [ OMP.ma_arr [ i ] ] , 2 ]
          := count_op [ OMP.op_arr [ i ] ] ;

          y := 1 ;

          repeat

            if ( young_arr [ y , 1 ] > 0 ) and ( young_arr [ y , 2 ] <= min_time )

            then

            begin

              if ( t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] < temp_test1 [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] ) then

              begin

                t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ]
                := t_test [ temp_source1_match [ young_arr [ y , 1 ] , 1 ] ] + 1 ;

              end ;

              if ( t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] < temp_hand1 [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] ) then

              begin

                t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ]
                := t_hand [ temp_source1_match [ young_arr [ y , 1 ] , 2 ] ] + 1 ;    //三种资源各加一

              end;

              if ( t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] < temp_accessory1 [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] ) then

              begin
                t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ]
                := t_acces [ temp_source1_match [ young_arr [ y , 1 ] , 3 ] ] + 1 ;
              end ;

              pr_set := pr_set - [ young_arr [ y , 1 ] ] ;

            end ;

            y := y + 1 ;

          until young_arr [ y , 1 ] = 0 ;

          t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ]
          := t_test [ temp_source1_match [ OMP.ma_arr [ i ] , 1 ] ] - 1 ;

          t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ]
          := t_hand [ temp_source1_match [ OMP.ma_arr [ i ] , 2 ] ] - 1 ;       //三种资源各减一

          t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ]
          := t_acces [ temp_source1_match [ OMP.ma_arr [ i ] , 3 ] ] - 1 ;

          if not ( OMP.ma_arr [ i ] in pr_set ) then                            //将机器放入正在运行机器集内

          begin

            pr_set := pr_set + [ OMP.ma_arr [ i ] ] ;

          end;

        end ;

      end ;

    end ;

    source_count [ i , 1 , 1 ] := t_test [ 1 ] ;

    source_count [ i , 1 , 2 ] := t_test [ 2 ] ;

    source_count [ i , 1 , 3 ] := t_test [ 3 ] ;

    source_count [ i , 2 , 1 ] := t_hand [ 1 ] ;

    source_count [ i , 2 , 2 ] := t_hand [ 2 ] ;

    source_count [ i , 2 , 3 ] := t_hand [ 3 ] ;

    source_count [ i , 3 , 1 ] := t_acces [ 1 ] ;

    source_count [ i , 3 , 2 ] := t_acces [ 2 ] ;

    source_count [ i , 3 , 3 ] := t_acces [ 3 ] ;

    source_count [ i , 3 , 4 ] := t_acces [ 4 ] ;

  end ;



//////////////////////////////////////////////////////////////////////////////////////

  fitness := m_final [ 1 ] ;

  for i := 2 to 36 do

  begin

    if fitness < m_final [ i ] then

    begin

      fitness := m_final [ i ] ;

    end;

  end;

  for i := 1 to max_i do

  begin

    o11_arr [ i ] := OMP.op_arr [ i ] ;

    m11_arr [ i ] := OMP.ma_arr [ i ] ;

    pri_arr [ i ] := temp_o_arr [ i ] ;

  end;

end ;

function swap ( var swap_arr : arr ) : arr ;                   //操作1：随机交换

var
  i , j , t : Integer ;

begin

  i := random ( max_i) + 1 ;

  repeat

    j := random ( max_i) + 1 ;

  until j <> i ;

  t := swap_arr [ i ] ;

  swap_arr [ i ] := swap_arr [ j ] ;

  swap_arr [ j ] := t ;

  Result := swap_arr ;

end ;

function adj_swap ( var adj_swap_arr : arr ) : arr ;              //操作2：前后交换

var
  i , t : Integer ;

begin

  i := random ( max_i - 1 ) + 1 ;

  t := adj_swap_arr [ i ] ;

  adj_swap_arr [ i ] := adj_swap_arr [ i + 1 ] ;

  adj_swap_arr [ i + 1 ] := t ;

  Result := adj_swap_arr ;

end ;

function pre_insert ( var insert_arr : arr ) : arr ;              //操作3：前插

var
  i , u , v , t : integer ;
  insert_arr1 : arr ;

begin

  for i:= 1 to max_i do

  begin

    insert_arr1 [ i ] := insert_arr [ i ] ;

  end ;

  u := random ( max_i) + 1 ;

  repeat

    v := random ( max_i) + 1 ;

  until v <> u ;

  if v < u then

  begin

    t := u ;

    u := v ;

    v := t ;

  end ;

  for i := 1 to u - 1 do

  begin

    insert_arr [ i ] := insert_arr1 [ i ] ;

  end ;

  for i := v + 1 to max_i do

  begin

    insert_arr [ i ] := insert_arr1 [ i ] ;

  end ;

  insert_arr [ u ] :=insert_arr1 [ v ] ;

  for i := u + 1 to v do

  begin

    insert_arr [ i ] := insert_arr1 [ i - 1 ] ;

  end ;

  result := insert_arr ;

end;

function post_insert ( var back_insert_arr : arr ) : arr ;         //操作4：后插

var
  i , u , v , t : integer ;
  insert_arr1 : arr ;

begin

  for i:= 1 to max_i do

  begin

    insert_arr1 [ i ] := back_insert_arr [ i ] ;

  end ;

  u := random ( max_i) + 1 ;

  repeat

    v := random ( max_i) + 1 ;

  until v <> u ;

  if v < u then

  begin

    t := u ;

    u := v ;

    v := t ;

  end ;

  for i := 1 to u - 1 do

  begin

    back_insert_arr [ i ] := insert_arr1 [ i ] ;

  end ;

  for i := v + 1 to max_i do

  begin

    back_insert_arr [ i ] := insert_arr1 [ i ] ;

  end ;

  back_insert_arr [ v ] :=insert_arr1 [ u ] ;

  for i := u to v - 1 do

  begin

    back_insert_arr [ i ] := insert_arr1 [ i + 1 ] ;

  end ;

  result := back_insert_arr ;

end;

function binding_swap ( var binding_swap_arr : arr ) : arr ;      //操作5：邻域交换

var
  i , t , p : Integer ;

begin

  i := random ( max_i - 3 ) + 1 ;

  t := binding_swap_arr [ i ] ;

  binding_swap_arr [ i ] := binding_swap_arr [ i + 3 ] ;

  binding_swap_arr [ i + 3 ] := t ;


  
  p := binding_swap_arr [ i + 1 ] ;

  binding_swap_arr [ i + 1 ] := binding_swap_arr [ i + 2 ] ;

  binding_swap_arr [ i + 2 ] := p ;

  Result := binding_swap_arr ;

end ;

function hside_swap ( var hside_swap_arr : arr ) : arr ;          //操作6：区块四点交换

var
  divi : Double ;
  d , i , j , t , k , v : Integer ;

begin

  divi := max_i/2 ;

  d := trunc ( divi ) ;

  i := Random ( d ) + 1 ;

  repeat

    j := random ( d ) + 1 ;

  until j <> i ;

  if j < i then

  begin

    t := i ;

    i := j ;

    j := t ;

  end;

  k := random ( d ) + 1 + d ;

  repeat

    v := random ( d ) + 1 + d ;

  until v <> k ;

  if v < k then

  begin

    t := k ;

    k := v ;

    v := t ;

  end;

  t := hside_swap_arr [ i ] ;

  hside_swap_arr [ i ] := hside_swap_arr [ k ] ;

  hside_swap_arr [ k ] := t ;

  t := hside_swap_arr [ j ] ;

  hside_swap_arr [ j ] := hside_swap_arr [ v ] ;

  hside_swap_arr [ v ] := t ;

  Result := hside_swap_arr ;

end ;

function inverse ( var inverse_arr : arr ) : arr ;              //操作7：区块转置

var
  divi : Double ;
  d , i , j , t , n : Integer ;
  core_arr : arr ;

begin

  i := random ( max_i) + 1 ;

  repeat

    j := random ( max_i) + 1 ;

  until j <> i ;

  if j < i then

  begin

    t := i ;

    i := j ;

    j := t ;

  end ;

  divi := ( j - i ) / 2 ;

  d := Trunc ( divi ) ;

  for n := 1 to ( j - i + 1 ) do

  begin

    core_arr [ n ] := inverse_arr [ i + n - 1 ] ;

  end;

  for n := 1 to d do

  begin

    t := core_arr [ n ] ;

    core_arr [ n ] := core_arr [ j - i + 1 - n + 1 ] ;

    core_arr [ j - i + 1 - n + 1 ] := t ;

  end;

  for n := 1 to ( j - i + 1 ) do

  begin

    inverse_arr [ i + n - 1 ] := core_arr [ n ]

  end;

  Result := inverse_arr ;

end ;

function segment_swap ( var segment_arr : arr ) : arr ;               //操作8：部分交换

var
  i , j , n , t : Integer ;

begin

  i := random ( max_i - 6 ) + 1 ;

  repeat

    j := random ( max_i - 6 ) + 1 ;

  until ( j <> i ) and (  ( j - 6 > i ) or ( i - 6 > j ) );

  for n := 1 to 6 do

  begin

    t := segment_arr [ i + n - 1 ] ;

    segment_arr [ i + n - 1 ] := segment_arr [ j + n - 1 ] ;

    segment_arr [ j + n - 1 ] := t ;

  end ;

  Result := segment_arr ;

end ;

procedure hl_individual_evaluated_and_bestindiv_generated ( hl_evaluated_pop : hls_indiv_contrast_pop ;
                                                                      test1 , hand1 , accessory1 : arr ;
                                                                              source1_match : two_arr4 ;
                                                                                setup1_time : two_arr1 ;
                                                                   M_mt_arr1 , P_mt_arr1 : three_arr ;
                                                                                 var eva_agg : double_arr ) ;
var
  i , j , k , y , w : Integer ;
  evaluate_op , record_op , o11_arr , m11_arr , pri_arr  : arr ;
  test_hls_record : hls_record_pop ;
  test_hls_store : array [ 1 .. 2*hls_popsize ] of hls_record ;
  fitness : Double ;
  record_va : hls_record ;

begin

  for i := 1 to hls_popsize do

  begin

    for j := 1 to max_i do

    begin

      evaluate_op [ j ] := pai0.pri_op [ j ] ;      //即将以初始化生成的最优低层个体作为评价标准进行高层种群的评价

    end;

    for j := 1 to hls_length do                     //每一个操作π0进行一次变换

    begin

      case hl_evaluated_pop [ i ].arr [ j ] of

        1: begin

             swap ( evaluate_op ) ;

             calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

             for w := 1 to hls_length do                             //变序12次完之后取最优

             begin

               test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ; // test_hls_record{高层操作序；评价完之后的工序排序；高层评价值}

             end ;

             for k := 1 to max_i do

             begin

               test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

             end;

             test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        2: begin

              adj_swap ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        3: begin

              pre_insert ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        4: begin

              post_insert ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        5: begin

              binding_swap ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        6: begin

              hside_swap ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        7: begin

              inverse ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        8: begin

              segment_swap ( evaluate_op ) ;

              calc_fitness ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end ;

      end ;

    end;

    for y := 2 to hls_length  do      //↓↓↓变换完12次之后找到其中最好的fit来计算这一条高层操作序列的评价值

    begin

      if test_hls_record [ 1 ].eva < test_hls_record [ y ].eva then

      begin

        record_va := test_hls_record [ y ] ;

        test_hls_record [ y ] := test_hls_record [ 1 ] ;

        test_hls_record [ 1 ] := record_va ;

      end;

    end;

    eva_agg [ i ] :=  test_hls_record [ 1 ].eva ;          //留作给高层个体附上评价值

    test_hls_store [ i ] := test_hls_record [ 1 ] ;        //留作选出这一代最优op

  end;

   //选出这一代最优序列↓↓↓↓↓↓↓↓↓↓↓↓↓

  for y := 2 to hls_popsize  do

  begin

    if test_hls_store [ 1 ].eva < test_hls_store [ y ].eva then

    begin

      record_va := test_hls_store [ y ] ;

      test_hls_store [ y ] := test_hls_store [ 1 ] ;

      test_hls_store [ 1 ] := record_va ;

    end;

  end;

  for i := 1 to max_i do

  begin

    record_op [ i ] := test_hls_store [ 1 ].op [ i ] ;

  end;

  calc_fitness ( record_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

  if litera_num =1 then

  begin

    for i := 1 to max_i do

    begin

      best_contrast_ind.pri_op [ i ] := pri_arr [ i ] ;

      best_contrast_ind.op [ i ] := o11_arr [ i ] ;

      best_contrast_ind.ma [ i ] := m11_arr [ i ] ;

    end;

    best_contrast_ind.fitness := fitness ;

  end;

  if litera_num > 1 then

  begin

    if best_contrast_ind.fitness > fitness then

    begin

      for i := 1 to max_i do

      begin

        best_contrast_ind.pri_op [ i ] := pri_arr [ i ] ;

        best_contrast_ind.op [ i ] := o11_arr [ i ] ;

        best_contrast_ind.ma [ i ] := m11_arr [ i ] ;

      end;

      best_contrast_ind.fitness := fitness ;

    end ;

  end ;
  
end ;

procedure hl_individual_evaluated_and_bestindiv_generated_compare ( hl_evaluated_pop : hls_indiv_contrast_pop ;
                                                                      test1 , hand1 , accessory1 : arr ;
                                                                              source1_match : two_arr4 ;
                                                                                setup1_time : two_arr1 ;
                                                                   M_mt_arr1 , P_mt_arr1 : three_arr ;
                                                                                 var eva_agg : double_arr ) ;
var
  i , j , k , y , w : Integer ;
  evaluate_op , record_op , o11_arr , m11_arr , pri_arr  : arr ;
  test_hls_record : hls_record_pop ;
  test_hls_store : array [ 1 .. 2*hls_popsize ] of hls_record ;
  fitness : Double ;
  record_va : hls_record ;

begin

  for i := 1 to hls_popsize do

  begin

    for j := 1 to max_i do

    begin

      evaluate_op [ j ] := pai0.pri_op [ j ] ;      //即将以初始化生成的最优低层个体作为评价标准进行高层种群的评价

    end;

    for j := 1 to hls_length do                     //每一个操作π0进行一次变换

    begin

      case hl_evaluated_pop [ i ].arr [ j ] of

        1: begin

             swap ( evaluate_op ) ;

             calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

             for w := 1 to hls_length do                             //变序12次完之后取最优

             begin

               test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ; // test_hls_record{高层操作序；评价完之后的工序排序；高层评价值}

             end ;

             for k := 1 to max_i do

             begin

               test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

             end;

             test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        2: begin

              adj_swap ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        3: begin

              pre_insert ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        4: begin

              post_insert ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        5: begin

              binding_swap ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        6: begin

              hside_swap ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        7: begin

              inverse ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end;

        8: begin

              segment_swap ( evaluate_op ) ;

              calc_fitness_compare ( evaluate_op , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

              for w := 1 to hls_length do

              begin

                test_hls_record [ j ].arr [ w ] := hl_evaluated_pop [ i ].arr [ w ] ;

              end ;

              for k := 1 to max_i do

              begin

                test_hls_record [ j ].op [ k ] := evaluate_op [ k ] ;

              end;

              test_hls_record [ j ].eva := pai0.fitness - fitness ;

           end ;

      end ;

    end;

    for y := 2 to hls_length  do      //↓↓↓变换完12次之后找到其中最好的fit来计算这一条高层操作序列的评价值

    begin

      if test_hls_record [ 1 ].eva < test_hls_record [ y ].eva then

      begin

        record_va := test_hls_record [ y ] ;

        test_hls_record [ y ] := test_hls_record [ 1 ] ;

        test_hls_record [ 1 ] := record_va ;

      end;

    end;

    eva_agg [ i ] :=  test_hls_record [ 1 ].eva ;          //留作给高层个体附上评价值

    test_hls_store [ i ] := test_hls_record [ 1 ] ;        //留作选出这一代最优op

  end;

   //选出这一代最优序列↓↓↓↓↓↓↓↓↓↓↓↓↓

  for y := 2 to hls_popsize  do

  begin

    if test_hls_store [ 1 ].eva < test_hls_store [ y ].eva then

    begin

      record_va := test_hls_store [ y ] ;

      test_hls_store [ y ] := test_hls_store [ 1 ] ;

      test_hls_store [ 1 ] := record_va ;

    end;

  end;

  for i := 1 to max_i do

  begin

    record_op [ i ] := test_hls_store [ 1 ].op [ i ] ;

  end;

  calc_fitness_compare ( record_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

  if litera_num =1 then

  begin

    for i := 1 to max_i do

    begin

      best_contrast_ind.pri_op [ i ] := pri_arr [ i ] ;

      best_contrast_ind.op [ i ] := o11_arr [ i ] ;

      best_contrast_ind.ma [ i ] := m11_arr [ i ] ;

    end;

    best_contrast_ind.fitness := fitness ;

  end;

  if litera_num > 1 then

  begin

    if best_contrast_ind.fitness > fitness then

    begin

      for i := 1 to max_i do

      begin

        best_contrast_ind.pri_op [ i ] := pri_arr [ i ] ;

        best_contrast_ind.op [ i ] := o11_arr [ i ] ;

        best_contrast_ind.ma [ i ] := m11_arr [ i ] ;

      end;

      best_contrast_ind.fitness := fitness ;

    end ;

  end ;

end ;

procedure Initialize_pop ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;
var
  y , i , j , k, pop_value , temper , w , q : Integer ;
  fitness : Double ;

  temp_o , trans_o_arr , temp_o_arr , copy_arr , o11_arr , m11_arr , pri_arr : arr ;
  eva_agg : double_arr ;
  t_set : set of 1 .. 255 ;
  flag : Boolean ;
  SBZM1 : SBZM ;

  hls_arr : array [ 1 .. hls_length ] of Integer ;
  count_hls , inibest_op : arr ;
  test_hls_record : hls_record_pop ;
  record_va : hls_record ;
  hl_record : hls_indiv ;

begin
  //////////以下是初始低层问题域种群、初始最好的个体生成////////////////////////
  randomize ;
  pop_value := initially_generated_individuals ;

  for i := 1 to max_i do

  begin

    copy_arr [ i ] := o1_arr [ i ] ;

  end ;

  for y := 1 to pop_value do

  begin

    temper := random ( max_i) + 1 ;

    temp_o [1] := temper ;

    t_set := [ ] ;

    t_set := t_set + [ temper]  ;

    for i := 2 to max_i do

    begin

      flag := True ;

      repeat

        temper := random ( max_i ) + 1 ;

        if not ( temper in t_set ) then

        begin

          temp_o [ i ] := temper ;

          t_set := t_set + [ temper ] ;

          flag := True ;

        end

        else

        begin

          flag := False ;

        end ;

      until flag = true ;

    end;

    for i := 1 to max_i do

    begin

      o1_arr [ i ] := copy_arr [ temp_o [ i ] ] ;

    end ;

    calc_fitness ( o1_arr , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

    SBZM_pop [ y ].fitness := fitness ;

    for w := 1 to max_i do

    begin

      SBZM_pop [ y ].op [ w ] := o11_arr [ w ] ;

      SBZM_pop [ y ].ma [ w ] := m11_arr [ w ] ;

      SBZM_pop [ y ].pri_op [ w ] := pri_arr [ w ] ;

    end ;

  end ;

/////////////////////////////////////////


  for i := 1 to pop_value - 1 do                         //initially_generated_individuals个随机生成的个体冒泡排序

  begin

    for j := i + 1 to pop_value do

    begin

      if sbzm_pop [ i ].fitness > sbzm_pop [ j ].fitness then

      begin

        sbzm1 := sbzm_pop [ i ] ;

        sbzm_pop [ i ] := sbzm_pop [ j ] ;

        sbzm_pop [ j ] := sbzm1 ;

      end;

    end ;

  end;

  for i := 1 to max_i do

  begin

    pai0.pri_op [ i ] := sbzm_pop [ 1 ].pri_op [ i ] ;

    pai0.op [ i ] := sbzm_pop [ 1 ].op [ i ] ;

    pai0.ma [ i ] := sbzm_pop [ 1 ].ma [ i ] ;

  end;

  pai0.fitness := sbzm_pop [ 1 ].fitness ;


 //////////以上是初始低层问题域////////////////////////////////////

 //////////以下是初始高层策略域////////////////////////////////////

  for i := 1 to hls_popsize do                                      //随机生成80个高层策略域个体的序列

  begin

    FillChar ( count_hls , sizeof ( arr ) , 0 ) ;

    temper := Random ( MC_width_and_hight ) + 1 ;

    hls_arr [ 1 ] := temper ;

    count_hls [ hls_arr [ 1 ] ] := count_hls [ hls_arr [ 1 ] ] + 1 ;

    for j := 2 to hls_length do

    begin

      flag := True ;

      repeat

        temper := random ( MC_width_and_hight ) + 1 ;

        if ( count_hls [ temper ] < 3 ) then

        begin

          hls_arr [ j ] := temper ;

          count_hls [ hls_arr [ j ] ] := count_hls [ hls_arr [ j ] ] + 1 ;

          flag := True ;

        end

        else
        begin

          flag := False ;

        end ;

      until flag = true ;

    end;

    for w := 1 to hls_length do

    begin

      hls_pop [ i ].arr [ w ] := hls_arr [ w ] ;

    end;
  end;

  /////////////下面开始对随机生成的高层个体进行评估值计算并取优/////////////////

  hl_individual_evaluated_and_bestindiv_generated ( hls_pop ,  test1 , hand1 , accessory1 ,
                                                                                    source1_match ,
                                                                                      setup1_time ,
                                                                         M_mt_arr1 , P_mt_arr1 , eva_agg ) ;

  for i := 1 to hls_popsize do

  begin

    hls_pop [ i ].evalue := eva_agg [ i ] ;

  end;

  for i := 1 to hls_popsize - 1 do                         //高层个体冒泡排序

  begin

    for j := i + 1 to hls_popsize do

    begin

      if hls_pop [ i ].evalue < hls_pop [ j ].evalue then

      begin

        hl_record := hls_pop [ i ] ;

        hls_pop [ i ] := hls_pop [ j ] ;

        hls_pop [ j ] := hl_record ;

      end;

    end ;

  end ;

end ;

procedure Initialize_pop_compare ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;
var
  y , i , j , k, pop_value , temper , w , q : Integer ;
  fitness : Double ;

  temp_o , trans_o_arr , temp_o_arr , copy_arr , o11_arr , m11_arr , pri_arr : arr ;
  eva_agg : double_arr ;
  t_set : set of 1 .. 255 ;
  flag : Boolean ;
  SBZM1 : SBZM ;

  hls_arr : array [ 1 .. hls_length ] of Integer ;
  count_hls , inibest_op : arr ;
  test_hls_record : hls_record_pop ;
  record_va : hls_record ;
  hl_record : hls_indiv ;

begin
  //////////以下是初始低层问题域种群、初始最好的个体生成////////////////////////
  randomize ;
  pop_value := initially_generated_individuals ;

  for i := 1 to max_i do

  begin

    copy_arr [ i ] := o1_arr [ i ] ;

  end ;

  for y := 1 to pop_value do

  begin

    temper := random ( max_i) + 1 ;

    temp_o [1] := temper ;

    t_set := [ ] ;

    t_set := t_set + [ temper]  ;

    for i := 2 to max_i do

    begin

      flag := True ;

      repeat

        temper := random ( max_i ) + 1 ;

        if not ( temper in t_set ) then

        begin

          temp_o [ i ] := temper ;

          t_set := t_set + [ temper ] ;

          flag := True ;

        end

        else

        begin

          flag := False ;

        end ;

      until flag = true ;

    end;

    for i := 1 to max_i do

    begin

      o1_arr [ i ] := copy_arr [ temp_o [ i ] ] ;

    end ;

    calc_fitness_compare ( o1_arr , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

    SBZM_pop [ y ].fitness := fitness ;

    for w := 1 to max_i do

    begin

      SBZM_pop [ y ].op [ w ] := o11_arr [ w ] ;

      SBZM_pop [ y ].ma [ w ] := m11_arr [ w ] ;

      SBZM_pop [ y ].pri_op [ w ] := pri_arr [ w ] ;

    end ;

  end ;

/////////////////////////////////////////


  for i := 1 to pop_value - 1 do                         //initially_generated_individuals个随机生成的个体冒泡排序

  begin

    for j := i + 1 to pop_value do

    begin

      if sbzm_pop [ i ].fitness > sbzm_pop [ j ].fitness then

      begin

        sbzm1 := sbzm_pop [ i ] ;

        sbzm_pop [ i ] := sbzm_pop [ j ] ;

        sbzm_pop [ j ] := sbzm1 ;

      end;

    end ;

  end;

  for i := 1 to max_i do

  begin

    pai0.pri_op [ i ] := sbzm_pop [ 1 ].pri_op [ i ] ;

    pai0.op [ i ] := sbzm_pop [ 1 ].op [ i ] ;

    pai0.ma [ i ] := sbzm_pop [ 1 ].ma [ i ] ;

  end;

  pai0.fitness := sbzm_pop [ 1 ].fitness ;


 //////////以上是初始低层问题域////////////////////////////////////

 //////////以下是初始高层策略域////////////////////////////////////

  for i := 1 to hls_popsize do                                      //随机生成80个高层策略域个体的序列

  begin

    FillChar ( count_hls , sizeof ( arr ) , 0 ) ;

    temper := Random ( MC_width_and_hight ) + 1 ;

    hls_arr [ 1 ] := temper ;

    count_hls [ hls_arr [ 1 ] ] := count_hls [ hls_arr [ 1 ] ] + 1 ;

    for j := 2 to hls_length do

    begin

      flag := True ;

      repeat

        temper := random ( MC_width_and_hight ) + 1 ;

        if ( count_hls [ temper ] < 3 ) then

        begin

          hls_arr [ j ] := temper ;

          count_hls [ hls_arr [ j ] ] := count_hls [ hls_arr [ j ] ] + 1 ;

          flag := True ;

        end

        else
        begin

          flag := False ;

        end ;

      until flag = true ;

    end;

    for w := 1 to hls_length do

    begin

      hls_pop [ i ].arr [ w ] := hls_arr [ w ] ;

    end;
  end;

  /////////////下面开始对随机生成的高层个体进行评估值计算并取优/////////////////

  hl_individual_evaluated_and_bestindiv_generated_compare ( hls_pop ,  test1 , hand1 , accessory1 ,
                                                                                    source1_match ,
                                                                                      setup1_time ,
                                                                         M_mt_arr1 , P_mt_arr1 , eva_agg ) ;

  for i := 1 to hls_popsize do

  begin

    hls_pop [ i ].evalue := eva_agg [ i ] ;

  end;

  for i := 1 to hls_popsize - 1 do                         //高层个体冒泡排序

  begin

    for j := i + 1 to hls_popsize do

    begin

      if hls_pop [ i ].evalue < hls_pop [ j ].evalue then

      begin

        hl_record := hls_pop [ i ] ;

        hls_pop [ i ] := hls_pop [ j ] ;

        hls_pop [ j ] := hl_record ;

      end;

    end ;

  end ;

end ;

procedure Initialize_PM_Matrix  ;           //初始化概率矩阵

var
  i , j , k   : Integer ;

begin

  fillchar ( Pro_model_MC , sizeof ( real_three_arr ) , 0 ) ;

  for k := 1 to MC_length do

  begin

    if ( k = 1 ) then

    begin

      for i := 1 to MC_width_and_hight do

      begin

        for j := 1 to MC_width_and_hight do

        begin

          Pro_model_MC [ k , i , j ] := 1 / ( MC_width_and_hight ) ;          //cube size（高层序列长度）更改时记得在type里面更改相应长度

        end;

      end;

    end

    else
    begin

      for i := 1 to MC_width_and_hight do

      begin

        for j := 1 to MC_width_and_hight do

        begin

          Pro_model_MC [ k , i , j ] := 1 / ( MC_width_and_hight * MC_width_and_hight ) ;

        end;

      end;

    end ;


  end;

end ;

procedure Creating_MC_Matrix ;           //生成相似块矩阵

var
  i , j , k : Integer ;


begin

  fillchar ( similar_blocks_MC , sizeof ( real_three_arr ) , 0 ) ;



  for i := 1 to hls_length do

  begin

    for j := 1 to MC_length do

    begin

      similar_blocks_MC [ j , hls_pop [ i ].arr [ j ] , hls_pop [ i ].arr [ j + 1 ] ] :=
      similar_blocks_MC [ j , hls_pop [ i ].arr [ j ] , hls_pop [ i ].arr [ j + 1 ] ] + 1 ;

    end;

  end;

end ;

procedure creating_1st_PM_Matrix ;           //生成第一代概率矩阵

var
  i , j , k  : Integer ;

begin

  for k := 1 to MC_length do

  begin

    if ( k = 1 ) then

    begin

      for i := 1 to MC_width_and_hight do

      begin

        for j := 1 to MC_width_and_hight do

        begin

          Pro_model_MC [ k , i , j ] := similar_blocks_MC [ k , i , j ] / ( hls_length ) ;

        end ;

      end ;

    end

    else
    begin

      for i := 1 to MC_width_and_hight do

      begin

        for j := 1 to MC_width_and_hight do

        begin

          Pro_model_MC [ k , i , j ] := ( similar_blocks_MC [ k , i , j ] + Pro_model_MC [ k , i , j ]  )
                                        / ( 1 + hls_length ) ;

        end ;

      end ;

    end ;

  end ;

end ;

procedure creating_subs_PM_Matrix ;

var
  i , j , k  : Integer ;

begin

  for k := 1 to MC_length do

  begin

    for i := 1 to MC_width_and_hight do

    begin

      for j := 1 to MC_width_and_hight do

      begin

        Pro_model_MC [ k , i , j ] := ( 1 - r ) * ( Pro_model_MC [ k , i , j ]) +
                                      r * ( similar_blocks_MC [ k , i , j ] ) / hls_length ;

      end;

    end;

  end;

end;

procedure Updating_HPOP_basedon_PMMC (  test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
  i , j , k , w : integer ;
  sum_p : real ;
  //point_value : real ;                 //概率指针
  sim_blc_psum_arr : psum_arr ;
  eva_agg : double_arr ;
  hl_record : hls_indiv ;

begin

  randomize ;

  for i := 1 to hls_popsize do             //开始生成新的高层种群（i个高层个体）

  begin

    for j := 1 to hls_length do            //j个位置

    begin

      if ( j = 1 ) then                    //第一个位置

      begin
        // 先生成point value和sim blc psum arr ，然后采样生成相对应位置的高层操作
        point_value := random ;

        sum_p := 0 ;

        for k := 1 to MC_width_and_hight do        //sump作为选定第一个高层操作序的概率数生成范围

        begin

          for w := 1 to MC_width_and_hight do

          begin

            sum_p := sum_p + Pro_model_MC [ j , k , w ] ;

          end;

        end;

        point_value := ( point_value ) * ( sum_p ) ;

        fillchar ( sim_blc_psum_arr , sizeof ( psum_arr ) , 0 ) ;

        for k := 1 to MC_width_and_hight do

        begin

          if ( k = 1 ) then

          begin

            for w := 1 to MC_width_and_hight do

            begin

              sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + Pro_model_MC [ j , k , w ] ;

            end;

          end

          else
          begin

            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k - 1 ] ;

            for w := 1 to MC_width_and_hight do

            begin

              sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + Pro_model_MC [ j , k , w ] ;

            end;

          end;

        end;

        //第一个位置从1~8操作选一个生成   ↓↓↓↓

        if ( point_value >= 0 ) and ( point_value <= sim_blc_psum_arr [ 1 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 1 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 1 ] ) and ( point_value <= sim_blc_psum_arr [ 2 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 2 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 2 ] ) and ( point_value <= sim_blc_psum_arr [ 3 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 3 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 3 ] ) and ( point_value <= sim_blc_psum_arr [ 4 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 4 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 4 ] ) and ( point_value <= sim_blc_psum_arr [ 5 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 5 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 5 ] ) and ( point_value <= sim_blc_psum_arr [ 6 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 6 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 6 ] ) and ( point_value <= sim_blc_psum_arr [ 7 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 7 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 7 ] ) and ( point_value <= sim_blc_psum_arr [ 8 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 8 ;

        end ;

      end

      else               //不是第一个位置
      begin
        // 先生成point value和sim blc psum arr ，然后采样生成相对应位置的高层操作
        point_value := random ;

        sum_p := 0 ;

        for w := 1 to MC_width_and_hight do        //sump作为选定高层操作序的概率数生成范围

        begin

          sum_p := sum_p + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , w ] ;

        end;

        point_value := ( point_value ) * ( sum_p ) ;

        fillchar ( sim_blc_psum_arr , sizeof ( psum_arr ) , 0 ) ;

        for k := 1 to MC_width_and_hight do

        begin

          if ( k = 1 ) then

          begin



              sim_blc_psum_arr [ k ] :=
              sim_blc_psum_arr [ k ] + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , k ] ;



          end

          else
          begin

            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k - 1 ] ;

            sim_blc_psum_arr [ k ] :=
            sim_blc_psum_arr [ k ] + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , k ] ;

          end;

        end;

        //第一个位置从1~8操作选一个生成   ↓↓↓↓

        if ( point_value >= 0 ) and ( point_value <= sim_blc_psum_arr [ 1 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 1 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 1 ] ) and ( point_value <= sim_blc_psum_arr [ 2 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 2 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 2 ] ) and ( point_value <= sim_blc_psum_arr [ 3 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 3 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 3 ] ) and ( point_value <= sim_blc_psum_arr [ 4 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 4 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 4 ] ) and ( point_value <= sim_blc_psum_arr [ 5 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 5 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 5 ] ) and ( point_value <= sim_blc_psum_arr [ 6 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 6 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 6 ] ) and ( point_value <= sim_blc_psum_arr [ 7 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 7 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 7 ] ) and ( point_value <= sim_blc_psum_arr [ 8 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 8 ;

        end ;

      end;

    end ;

  end ;

  /////////////下面开始对高层个体进行评估值计算并取优/////////////////

  hl_individual_evaluated_and_bestindiv_generated ( hls_contrast_pop , test1 , hand1 , accessory1 ,
                                                                                    source1_match ,
                                                                                      setup1_time ,
                                                                         M_mt_arr1 , P_mt_arr1 , eva_agg ) ;

  for i := 1 to hls_popsize do                                 //  赋上评估值

  begin

    hls_contrast_pop [ i ].evalue := eva_agg [ i ] ;

  end;

  for i := hls_popsize + 1 to hls_popsize + hls_length do     //加上目前最优的高层

  begin

    hls_contrast_pop [ i ] := hls_pop [ i - hls_popsize ] ;

  end;

  for i := 1 to hls_popsize + hls_length - 1 do                         //总的高层个体冒泡排序

  begin

    for j := i + 1 to hls_popsize + hls_length do

    begin

      if hls_contrast_pop [ i ].evalue < hls_contrast_pop [ j ].evalue then

      begin

        hl_record := hls_contrast_pop [ i ] ;

        hls_contrast_pop [ i ] := hls_contrast_pop [ j ] ;

        hls_contrast_pop [ j ] := hl_record ;

      end;

    end ;

  end ;

  for i := 1 to hls_length do                                     //取前x个保存用作生成相似块矩阵

  begin

    hls_pop [ i ] := hls_contrast_pop [ i ] ;

  end;

end ;

procedure Updating_HPOP_basedon_PMMC_compare (  test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
  i , j , k , w : integer ;
  sum_p : real ;
  //point_value : real ;                 //概率指针
  sim_blc_psum_arr : psum_arr ;
  eva_agg : double_arr ;
  hl_record : hls_indiv ;

begin

  randomize ;

  for i := 1 to hls_popsize do             //开始生成新的高层种群（i个高层个体）

  begin

    for j := 1 to hls_length do            //j个位置

    begin

      if ( j = 1 ) then                    //第一个位置

      begin
        // 先生成point value和sim blc psum arr ，然后采样生成相对应位置的高层操作
        point_value := random ;

        sum_p := 0 ;

        for k := 1 to MC_width_and_hight do        //sump作为选定第一个高层操作序的概率数生成范围

        begin

          for w := 1 to MC_width_and_hight do

          begin

            sum_p := sum_p + Pro_model_MC [ j , k , w ] ;

          end;

        end;

        point_value := ( point_value ) * ( sum_p ) ;

        fillchar ( sim_blc_psum_arr , sizeof ( psum_arr ) , 0 ) ;

        for k := 1 to MC_width_and_hight do

        begin

          if ( k = 1 ) then

          begin

            for w := 1 to MC_width_and_hight do

            begin

              sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + Pro_model_MC [ j , k , w ] ;

            end;

          end

          else
          begin

            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k - 1 ] ;

            for w := 1 to MC_width_and_hight do

            begin

              sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + Pro_model_MC [ j , k , w ] ;

            end;

          end;

        end;

        //第一个位置从1~8操作选一个生成   ↓↓↓↓

        if ( point_value >= 0 ) and ( point_value <= sim_blc_psum_arr [ 1 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 1 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 1 ] ) and ( point_value <= sim_blc_psum_arr [ 2 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 2 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 2 ] ) and ( point_value <= sim_blc_psum_arr [ 3 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 3 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 3 ] ) and ( point_value <= sim_blc_psum_arr [ 4 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 4 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 4 ] ) and ( point_value <= sim_blc_psum_arr [ 5 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 5 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 5 ] ) and ( point_value <= sim_blc_psum_arr [ 6 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 6 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 6 ] ) and ( point_value <= sim_blc_psum_arr [ 7 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 7 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 7 ] ) and ( point_value <= sim_blc_psum_arr [ 8 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 8 ;

        end ;

      end

      else               //不是第一个位置
      begin
        // 先生成point value和sim blc psum arr ，然后采样生成相对应位置的高层操作
        point_value := random ;

        sum_p := 0 ;

        for w := 1 to MC_width_and_hight do        //sump作为选定高层操作序的概率数生成范围

        begin

          sum_p := sum_p + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , w ] ;

        end;

        point_value := ( point_value ) * ( sum_p ) ;

        fillchar ( sim_blc_psum_arr , sizeof ( psum_arr ) , 0 ) ;

        for k := 1 to MC_width_and_hight do

        begin

          if ( k = 1 ) then

          begin



              sim_blc_psum_arr [ k ] :=
              sim_blc_psum_arr [ k ] + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , k ] ;



          end

          else
          begin

            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k - 1 ] ;

            sim_blc_psum_arr [ k ] :=
            sim_blc_psum_arr [ k ] + Pro_model_MC [ j - 1 , hls_contrast_pop [ i ].arr [ j - 1 ] , k ] ;

          end;

        end;

        //第一个位置从1~8操作选一个生成   ↓↓↓↓

        if ( point_value >= 0 ) and ( point_value <= sim_blc_psum_arr [ 1 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 1 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 1 ] ) and ( point_value <= sim_blc_psum_arr [ 2 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 2 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 2 ] ) and ( point_value <= sim_blc_psum_arr [ 3 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 3 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 3 ] ) and ( point_value <= sim_blc_psum_arr [ 4 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 4 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 4 ] ) and ( point_value <= sim_blc_psum_arr [ 5 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 5 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 5 ] ) and ( point_value <= sim_blc_psum_arr [ 6 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 6 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 6 ] ) and ( point_value <= sim_blc_psum_arr [ 7 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 7 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 7 ] ) and ( point_value <= sim_blc_psum_arr [ 8 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 8 ;

        end ;

      end;

    end ;

  end ;

  /////////////下面开始对高层个体进行评估值计算并取优/////////////////

  hl_individual_evaluated_and_bestindiv_generated_compare ( hls_contrast_pop , test1 , hand1 , accessory1 ,
                                                                                    source1_match ,
                                                                                      setup1_time ,
                                                                         M_mt_arr1 , P_mt_arr1 , eva_agg ) ;

  for i := 1 to hls_popsize do                                 //  赋上评估值

  begin

    hls_contrast_pop [ i ].evalue := eva_agg [ i ] ;

  end;

  for i := hls_popsize + 1 to hls_popsize + hls_length do     //加上目前最优的高层

  begin

    hls_contrast_pop [ i ] := hls_pop [ i - hls_popsize ] ;

  end;

  for i := 1 to hls_popsize + hls_length - 1 do                         //总的高层个体冒泡排序

  begin

    for j := i + 1 to hls_popsize + hls_length do

    begin

      if hls_contrast_pop [ i ].evalue < hls_contrast_pop [ j ].evalue then

      begin

        hl_record := hls_contrast_pop [ i ] ;

        hls_contrast_pop [ i ] := hls_contrast_pop [ j ] ;

        hls_contrast_pop [ j ] := hl_record ;

      end;

    end ;

  end ;

  for i := 1 to hls_length do                                     //取前x个保存用作生成相似块矩阵

  begin

    hls_pop [ i ] := hls_contrast_pop [ i ] ;

  end;

end ;


procedure main ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
 best_o : arr ;
 i , q : integer ;
 o11_arr , m11_arr , pri_arr : arr ;
 fitness , x , avg : Double ;
 gantt : gant_tu ;

begin

  x := 0 ;

  q := 1 ;

  repeat

    fillchar ( best_fitness_arr , sizeof ( arr_max ) , 0 ) ;

    time1 := gettickcount ;

    litera_num := 1 ;

    Initialize_pop ( o1_arr , test1 , hand1 , accessory1 ,        //工序；三种资源数目
                                           source1_match ,        //机器相匹配资源种类
                                             setup1_time ,        //机器间设置时间
                                 M_mt_arr1 , P_mt_arr1 ) ;        //机器和加工时间矩阵

    Initialize_PM_Matrix ;                       //初始化三维概率矩阵

    repeat

      if ( litera_num = 1 ) then

      begin

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      if ( litera_num = 2 ) then
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_1st_PM_Matrix ;                     //生成第一代三维概率矩阵

        Updating_HPOP_basedon_PMMC (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_subs_PM_Matrix ;                     //生成后续三维概率矩阵

        Updating_HPOP_basedon_PMMC (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end;

      time3 := gettickcount - time1 ;

    until time3 > max_i*36 ;

    for i := 1 to max_i do

    begin

      best_o [ i ] := best_contrast_ind.pri_op [ i ] ;

    end;

    calc_fit_and_data_stored ( best_o , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

    canshu_arr [ q ] := fitness ;

    x := x + fitness ;

    for i := 1 to max_i do

    begin

      gantt.start [ i ] := start_time [ i ] ;

      gantt.dura [ i ] := dura_time [ i ] ;

      gantt.macInfo [ i ] := m11_arr [ i ] ;

      gantt.jobId [ i ] := o11_arr [ i ] ;

      gantt.oper [ i ] := count_xu [ i ] ;

    end;

    q := q + 1 ;

  until ( q>1 ) ;



end ;

procedure NFOA ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
 best_o : arr ;
 i , q : integer ;
 o11_arr , m11_arr , pri_arr : arr ;
 fitness , x , avg : Double ;
 gantt : gant_tu ;

begin

  x := 0 ;

  q := 1 ;

  repeat

    fillchar ( best_fitness_arr , sizeof ( arr_max ) , 0 ) ;

    time1 := gettickcount ;

    litera_num := 1 ;

    Initialize_pop_compare ( o1_arr , test1 , hand1 , accessory1 ,        //工序；三种资源数目
                                           source1_match ,        //机器相匹配资源种类
                                             setup1_time ,        //机器间设置时间
                                 M_mt_arr1 , P_mt_arr1 ) ;        //机器和加工时间矩阵

    Initialize_PM_Matrix ;                       //初始化三维概率矩阵

    repeat

      if ( litera_num = 1 ) then

      begin

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      if ( litera_num = 2 ) then
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_1st_PM_Matrix ;                     //生成第一代三维概率矩阵

        Updating_HPOP_basedon_PMMC_compare (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_subs_PM_Matrix ;                     //生成后续三维概率矩阵

        Updating_HPOP_basedon_PMMC_compare (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end;

      time3 := gettickcount - time1 ;

    until time3 > max_i*36 ;

    {for i := 1 to max_i do

    begin

      best_o [ i ] := best_contrast_ind.pri_op [ i ] ;

    end;

    calc_fit_and_data_stored ( best_o , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

    canshu_arr [ q ] := fitness ;

    x := x + fitness ;

    for i := 1 to max_i do

    begin

      gantt.start [ i ] := start_time [ i ] ;

      gantt.dura [ i ] := dura_time [ i ] ;

      gantt.macInfo [ i ] := m11_arr [ i ] ;

      gantt.jobId [ i ] := o11_arr [ i ] ;

      gantt.oper [ i ] := count_xu [ i ] ;

    end; }

    canshu_arr [ q ] := best_contrast_ind.fitness ;

    x := x + best_contrast_ind.fitness ;

    q := q + 1 ;

  until ( q>25 ) ;



end ;

procedure HHMEDA_compare1 ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
 best_o : arr ;
 i , q : integer ;
 o11_arr , m11_arr , pri_arr : arr ;
 fitness , x , avg : Double ;
 gantt : gant_tu ;

begin

  x := 0 ;

  q := 1 ;

  repeat

    fillchar ( best_fitness_arr , sizeof ( arr_max ) , 0 ) ;

    time1 := gettickcount ;

    litera_num := 1 ;

    Initialize_pop_compare ( o1_arr , test1 , hand1 , accessory1 ,        //工序；三种资源数目
                                           source1_match ,        //机器相匹配资源种类
                                             setup1_time ,        //机器间设置时间
                                 M_mt_arr1 , P_mt_arr1 ) ;        //机器和加工时间矩阵

    Initialize_PM_Matrix ;                       //初始化三维概率矩阵

    repeat

      if ( litera_num = 1 ) then

      begin

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      if ( litera_num = 2 ) then
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_1st_PM_Matrix ;                     //生成第一代三维概率矩阵

        Updating_HPOP_basedon_PMMC_compare (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end

      else
      begin

        Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

        creating_subs_PM_Matrix ;                     //生成后续三维概率矩阵

        Updating_HPOP_basedon_PMMC_compare (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                               source1_match ,
                                                 setup1_time ,
                                     M_mt_arr1 , P_mt_arr1 ) ;

        if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

        begin

          best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

        end ;

        best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

        litera_num := litera_num + 1 ;

      end;

      time3 := gettickcount - time1 ;

    until time3 > max_i*36 ;

    {for i := 1 to max_i do

    begin

      best_o [ i ] := best_contrast_ind.pri_op [ i ] ;

    end;

    calc_fit_and_data_stored ( best_o , test1 , hand1 , accessory1 ,
                                         source1_match ,
                                           setup1_time ,
                                 M_mt_arr1 , P_mt_arr1 ,
                                    o11_arr , m11_arr , pri_arr , fitness ) ;

    canshu_arr [ q ] := fitness ;

    x := x + fitness ;

    for i := 1 to max_i do

    begin

      gantt.start [ i ] := start_time [ i ] ;

      gantt.dura [ i ] := dura_time [ i ] ;

      gantt.macInfo [ i ] := m11_arr [ i ] ;

      gantt.jobId [ i ] := o11_arr [ i ] ;

      gantt.oper [ i ] := count_xu [ i ] ;

    end; }

    canshu_arr [ q ] := best_contrast_ind.fitness ;

    x := x + best_contrast_ind.fitness ;

    q := q + 1 ;

  until ( q>25 ) ;



end ;

procedure Initialize_EDA_MC  ;           //初始化EDA采样矩阵

var
  i , j : Integer ;

begin

  fillchar ( EDA_MC , sizeof ( real_two_arr ) , 0 ) ;

  for i := 1 to hls_length do

  begin

    for j := 1 to MC_width_and_hight do

    begin

      EDA_MC [ i , j ] := 1 / ( MC_width_and_hight * MC_width_and_hight ) ;

    end;

  end;

end ;

procedure Updeating_EDA_MC  ;           ////更新EDA采样矩阵

var
  i , j : Integer ;

begin

  Initialize_EDA_MC ;

  for i := 1 to hls_length do

  begin

    for j := 1 to hls_length do

    begin

      EDA_MC [ i , hls_pop [ i ].arr [ j ] ] :=
      EDA_MC [ i , hls_pop [ i ].arr [ j ] ] + 1 ;

    end;

  end;

end ;

procedure Updating_HPOP_basedon_EDAMC (  test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
  i , j , k , w : integer ;
  sum_p : real ;
  //point_value : real ;                 //概率指针
  sim_blc_psum_arr : psum_arr ;
  eva_agg : double_arr ;
  hl_record : hls_indiv ;

begin

  randomize ;

  for i := 1 to hls_popsize do             //开始生成新的高层种群（i个高层个体）

  begin

    for j := 1 to hls_length do            //j个位置

    begin


        // 先生成point value和sim blc psum arr ，然后采样生成相对应位置的高层操作
        point_value := random ;

        sum_p := 0 ;

        for k := 1 to MC_width_and_hight do        //sump作为选定第一个高层操作序的概率数生成范围

        begin



            sum_p := sum_p + EDA_MC [ j , k ] ;



        end;

        point_value := ( point_value ) * ( sum_p ) ;

        fillchar ( sim_blc_psum_arr , sizeof ( psum_arr ) , 0 ) ;

        for k := 1 to MC_width_and_hight do

        begin

          if ( k = 1 ) then

          begin



              sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + EDA_MC [ j , k ] ;



          end

          else
          begin

            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k - 1 ] ;



            sim_blc_psum_arr [ k ] := sim_blc_psum_arr [ k ] + EDA_MC [ j , k ] ;



          end;

        end;

        //第一个位置从1~8操作选一个生成   ↓↓↓↓

        if ( point_value >= 0 ) and ( point_value <= sim_blc_psum_arr [ 1 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 1 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 1 ] ) and ( point_value <= sim_blc_psum_arr [ 2 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 2 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 2 ] ) and ( point_value <= sim_blc_psum_arr [ 3 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 3 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 3 ] ) and ( point_value <= sim_blc_psum_arr [ 4 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 4 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 4 ] ) and ( point_value <= sim_blc_psum_arr [ 5 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 5 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 5 ] ) and ( point_value <= sim_blc_psum_arr [ 6 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 6 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 6 ] ) and ( point_value <= sim_blc_psum_arr [ 7 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 7 ;

        end

        else
        if ( point_value > sim_blc_psum_arr [ 7 ] ) and ( point_value <= sim_blc_psum_arr [ 8 ] ) then

        begin

          hls_contrast_pop [ i ].arr [ j ] := 8 ;

        end ;





    end ;

  end ;

  /////////////下面开始对高层个体进行评估值计算并取优/////////////////

  hl_individual_evaluated_and_bestindiv_generated ( hls_contrast_pop , test1 , hand1 , accessory1 ,
                                                                                    source1_match ,
                                                                                      setup1_time ,
                                                                         M_mt_arr1 , P_mt_arr1 , eva_agg ) ;

  for i := 1 to hls_popsize do                                 //  赋上评估值

  begin

    hls_contrast_pop [ i ].evalue := eva_agg [ i ] ;

  end;

  for i := hls_popsize + 1 to hls_popsize + hls_length do     //加上目前最优的高层

  begin

    hls_contrast_pop [ i ] := hls_pop [ i - hls_popsize ] ;

  end;

  for i := 1 to hls_popsize + hls_length - 1 do                         //总的高层个体冒泡排序

  begin

    for j := i + 1 to hls_popsize + hls_length do

    begin

      if hls_contrast_pop [ i ].evalue < hls_contrast_pop [ j ].evalue then

      begin

        hl_record := hls_contrast_pop [ i ] ;

        hls_contrast_pop [ i ] := hls_contrast_pop [ j ] ;

        hls_contrast_pop [ j ] := hl_record ;

      end;

    end ;

  end ;

  for i := 1 to hls_length do                                     //取前x个保存用作生成相似块矩阵

  begin

    hls_pop [ i ] := hls_contrast_pop [ i ] ;

  end;

end ;

procedure HHEDA ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
 best_o : arr ;
 i : integer ;
 o11_arr , m11_arr , pri_arr : arr ;
 fitness : Double ;
 gantt : gant_tu ;

begin

  litera_num := 1 ;

  Initialize_pop ( o1_arr , test1 , hand1 , accessory1 ,        //工序；三种资源数目
                                         source1_match ,        //机器相匹配资源种类
                                           setup1_time ,        //机器间设置时间
                               M_mt_arr1 , P_mt_arr1 ) ;        //机器和加工时间矩阵

  repeat

    if ( litera_num = 1 ) then

    begin

      best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

      litera_num := litera_num + 1 ;

    end

    else
    begin

      Creating_MC_Matrix ;                         //根据上一代种群生成三维相似块矩阵

      creating_subs_PM_Matrix ;                     //生成后续三维概率矩阵

      Updating_HPOP_basedon_EDAMC (  test1 , hand1 , accessory1 ,        //根据概率矩阵生成下一代高层种群
                                             source1_match ,
                                               setup1_time ,
                                   M_mt_arr1 , P_mt_arr1 ) ;

      if ( best_contrast_ind.fitness > best_fitness_arr [ litera_num - 1 ] ) then

      begin

        best_contrast_ind.fitness := best_fitness_arr [ litera_num - 1 ] ;

      end ;

      best_fitness_arr [ litera_num ] := best_contrast_ind.fitness ;

      litera_num := litera_num + 1 ;

    end;

  until litera_num = 99 ;

  for i := 1 to max_i do

  begin

    best_o [ i ] := best_contrast_ind.pri_op [ i ] ;

  end;

  calc_fit_and_data_stored ( best_o , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

  for i := 1 to max_i do

  begin

    gantt.start [ i ] := start_time [ i ] ;

    gantt.dura [ i ] := dura_time [ i ] ;

    gantt.macInfo [ i ] := m11_arr [ i ] ;

    gantt.jobId [ i ] := o11_arr [ i ] ;

    gantt.oper [ i ] := count_xu [ i ] ;

  end;

end ;

procedure QHH_state_choose ;

begin

  if ( QHH_pt >= 0 ) and ( QHH_pt < 0.85 ) then

  begin

    QHH_state := 1 ;

  end;

  if ( QHH_pt >= 0.85 ) and ( QHH_pt < 1 ) then

  begin

    QHH_state := 2 ;

  end;

  if ( QHH_pt >= 1 ) then

  begin

    QHH_state := 3 ;

  end;

end;

procedure QHH_action_choose ;

var
  k , St_a : Double ;
  i  : Integer ;

begin

  Randomize ;

  k := Random ;



  if k < QHH_greedy then

  begin

    QHH_action := Random ( 8 ) + 1 ;

  end
  else
  begin

    St_a := QHH_qtable [ QHH_state , 1 ] ;

    QHH_action := 1 ;

    for i := 2 to QHH_OP do

    begin

      if QHH_qtable [ QHH_state , i ] > St_a then

      begin

        St_a := QHH_qtable [ QHH_state , i ] ;

        QHH_action := i ;

      end;

    end;

  end;

end;

procedure QHH_perform_action_to_paib ( test1 , hand1 , accessory1 : arr ;
                                        source1_match : two_arr4 ;
                                          setup1_time : two_arr1 ;
                             M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
  i : Integer ;
  perdormed_op , o11_arr , m11_arr , pri_arr : arr ;
  fitness : Double ;

begin

  for i := 1 to max_i do

  begin

    perdormed_op [ i ] := pai_b.pri_op [ i ] ;

  end;



  begin

    case QHH_action of

      1: begin

           swap ( perdormed_op ) ;

           calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                     source1_match ,
                                       setup1_time ,
                             M_mt_arr1 , P_mt_arr1 ,
                                o11_arr , m11_arr , pri_arr , fitness ) ;

           if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      2: begin

            adj_swap ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

            if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      3: begin

            pre_insert ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

             if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      4: begin

            post_insert ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

             if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      5: begin

            binding_swap ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

            if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      6: begin

            hside_swap ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

            if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      7: begin

            inverse ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness) ;

            if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end;

      8: begin

            segment_swap ( perdormed_op ) ;

            calc_fitness ( perdormed_op , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

            if fitness < pai_b.fitness then

            begin

              for i := 1 to max_i do

              begin

                pai_b.pri_op [ i ] := pri_arr [ i ] ;

                pai_b.op [ i ] := o11_arr [ i ] ;

                pai_b.ma [ i ] := m11_arr  [ i ] ;

              end;

              pai_b.fitness := fitness ;

            end;

         end ;

    end ;

  end;



end;

procedure QHH_rein_sig_obtian ;

var
  u : Double ;

begin

  Randomize ;

  u := Random ;

  if ( QHH_greedy >= 0 ) and ( QHH_greedy < 0.85 ) then

  begin

    if ( u < QHH_greedy ) then

    begin

      QHH_rein_sig := 2 ;

    end

    else
    if ( u > QHH_greedy ) then

    begin

      QHH_rein_sig := 1 ;

    end

    else
    if ( u = QHH_greedy ) then

    begin

      QHH_rein_sig := 0 ;

    end;

  end;

  if ( QHH_greedy >= 0.85 ) and ( QHH_greedy < 1 ) then

  begin

    if ( u < QHH_greedy ) then

    begin

      QHH_rein_sig := 1 ;

    end

    else
    if ( u > QHH_greedy ) then

    begin

      QHH_rein_sig := 2 ;

    end

    else
    if ( u = QHH_greedy ) then

    begin

      QHH_rein_sig := 0 ;

    end;

  end;

end;

procedure QHH_argmax_v_calculate ;

var
  i : Integer ;

begin

  QHH_argmax_v := QHH_qtable [ QHH_state , 1 ] ;

  for i := 2 to QHH_OP do

  begin

    if QHH_qtable [ QHH_state , i ] > QHH_argmax_v then

    begin

      QHH_argmax_v := QHH_qtable [ QHH_state , i ] ;

    end;

  end;

end;

procedure QHH_qtable_updeated ;

var
  i , j : Integer ;

begin

  for i := 1 to 3 do

  begin

    for j := 1 to QHH_OP do

    begin

      QHH_qtable [ i , j ] := QHH_qtable [ i , j ] +
                              QHH_learn_rate*( QHH_rein_sig + QHH_discount*QHH_argmax_v - QHH_qtable [ i , j ] ) ;

    end;

  end;

end;

procedure QHH ( o1_arr , test1 , hand1 , accessory1 : arr ;
                                  source1_match : two_arr4 ;
                                    setup1_time : two_arr1 ;
                       M_mt_arr1 , P_mt_arr1 : three_arr ) ;

var
 best_o : arr ;
 i , j : integer ;
 o11_arr , m11_arr , pri_arr : arr ;
 fitness : Double ;
 gantt : gant_tu ;



begin

  fillchar ( best_fitness_arr , sizeof ( arr_max ) , 0 ) ;

  repeat
  QHH_G := 2200 ;

  QHH_EP := 8 ;

  QHH_discount := 0.4 ;

  Initialize_pop ( o1_arr , test1 , hand1 , accessory1 ,        //工序；三种资源数目
                                         source1_match ,        //机器相匹配资源种类
                                           setup1_time ,        //机器间设置时间
                               M_mt_arr1 , P_mt_arr1 ) ;        //机器和加工时间矩阵

  for i := 1 to max_i do                            //初始化πb

  begin

    pai_b.pri_op [ i ] := pai0.pri_op [ i ] ;

    pai_b.op [ i ] := pai0.op [ i ] ;

    pai_b.ma [ i ] := pai0.ma [ i ] ;

  end;

  pai_b.fitness := pai0.fitness ;

  QHH_pt := 1  ;                                   //初始化pt

  QHH_state_choose ;                              //初始化状态                 {调用函数}

  for i := 1 to 3 do                               //初始化q表

  begin

    for j := 1 to QHH_OP do

    begin

      QHH_qtable [ i , j ] := 0 ;

    end;

  end;

  QHH_Gcur := 1 ;

  repeat

    QHH_E := 1 ;

    QHH_greedy := 1 - ( QHH_gcur / QHH_G ) ;

    QHH_action_choose ;                                                         {调用函数}

    QHH_C_init := pai_b.fitness ;

    repeat

      QHH_perform_action_to_paib ( test1 , hand1 , accessory1 ,
                                  source1_match ,
                                    setup1_time ,
                       M_mt_arr1 , P_mt_arr1  );                                              {调用函数}

      QHH_E := QHH_E + 1 ;

    until QHH_E = QHH_EP ;

    QHH_C_EP := pai_b.fitness ;

    QHH_pt := QHH_C_init/QHH_C_EP ;

    QHH_state_choose ;                                                          {调用函数}

    QHH_rein_sig_obtian ;                                                       {调用函数 }

    QHH_argmax_v_calculate ;                                                    {调用函数 }

    QHH_learn_rate := 1 - ( ( 0.9*QHH_Gcur )/QHH_G ) ;

    QHH_qtable_updeated ;                                                       {调用函数 }

    best_fitness_arr [ QHH_Gcur ] := pai_b.fitness ;

    QHH_Gcur := QHH_Gcur + 1 ;

  until QHH_Gcur = QHH_G ;

  for i := 1 to max_i do

  begin

    best_o [ i ] := pai_b.pri_op [ i ] ;

  end;

  calc_fit_and_data_stored ( best_o , test1 , hand1 , accessory1 ,
                                       source1_match ,
                                         setup1_time ,
                               M_mt_arr1 , P_mt_arr1 ,
                                  o11_arr , m11_arr , pri_arr , fitness ) ;

  for i := 1 to max_i do

  begin

    gantt.start [ i ] := start_time [ i ] ;

    gantt.dura [ i ] := dura_time [ i ] ;

    gantt.macInfo [ i ] := m11_arr [ i ] ;

    gantt.jobId [ i ] := o11_arr [ i ] ;

    gantt.oper [ i ] := count_xu [ i ] ;

  end;

  until fitness = 180 ;
end ;



end.
