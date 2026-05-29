/*                                                                             
               SAS DATA DEFINITION STATEMENTS FOR ICPSR 9752                   
                    SURVEY OF CONSUMER FINANCES, 1977                          
                          SEPTEMBER 1992 VERSION                               
                                                                               
   The following SAS setup sections appear in this file for the LRECL          
   version of this data collection.  These sections are listed below:          
                                                                               
   PROC FORMAT:  creates user-defined formats for the variables.  Formats      
   replace original value codes with the value code descriptions.  Not all     
   variables necessarily have user-defined formats.                            
                                                                               
   DATA begins a SAS data step and names an output SAS dataset.                
                                                                               
   INFILE identifies the input file to be read with the input statement.       
   Users must replace the "file-specification" with host computer-specific     
   input file specifications.                                                  
                                                                               
   INPUT contains the SAS statements which assign the variable names and       
   specify the beginning and ending column locations in the LRECL data         
   file for each variable.                                                     
                                                                               
   LABEL assigns variable labels for all variables in the data file.           
   Variable labels and variable names may be identical for some data           
   files.                                                                      
                                                                               
   MISSING VALUE RECODE sets user-defined numeric missing values to            
   missing as interpreted by the SAS system.  Only variables with              
   user-defined missing values are included in these statements.               
                                                                               
   Users may combine and modify these sections or parts of these sections      
   to suit their specific needs.  Users will also need to change the           
   file-specification in the INFILE statement to an appropriate filename       
   for their system.  Please note that the MISSING VALUE RECODE section        
   has been commented out (i.e., '/*').  To include the MISSING VALUE          
   RECODE section in the final SAS setup, remove the comment indicators        
   from the section.                                                           
                                                                               
**************************************************************************** */
                                                                               
* SAS DATA, INFILE, INPUT STATEMENTS;                                          
                                                                               
DATA;                                                                          
INFILE "file-specification" LRECL=1294;                                        
INPUT                                                                          
   V1 1-4                   V3 5-5                   V4 6-7                    
   V5 8-9                   V6 10-10                 V7 11-11                  
   V8 12-12                 V9 13-13                 V10 14-14                 
   V11 15-15                V12 16-16                V13 17-17                 
   V14 18-18                V15 19-19                V16 20-20                 
   V17 21-21                V18 22-22                V19 23-23                 
   V20 24-24                V21 25-26                V22 27-28                 
   V23 29-30                V24 31-32                V25 33-34                 
   V26 35-35                V27 36-36                V28 37-37                 
   V29 38-38                V30 39-39                V31 40-40                 
   V32 41-41                V33 42-42                V34 43-43                 
   V35 44-44                V36 45-45                V37 46-46                 
   V38 47-47                V39 48-48                V40 49-49                 
   V41 50-50                V42 51-51                V43 52-52                 
   V44 53-53                V45 54-54                V46 55-55                 
   V47 56-56                V48 57-57                V49 58-58                 
   V50 59-59                V51 60-61                V52 62-63                 
   V53 64-65                V54 66-66                V55 67-67                 
   V56 68-68                V57 69-69                V58 70-70                 
   V59 71-72                V60 73-74                V61 75-75                 
   V62 76-76                V63 77-77                V64 78-78                 
   V65 79-79                V66 80-80                V67 81-81                 
   V68 82-82                V69 83-83                V70 84-85                 
   V71 86-87                V72 88-91                V73 92-92                 
   V74 93-94                V75 95-95                V76 96-96                 
   V77 97-97                V78 98-98                V79 99-99                 
   V80 100-100              V81 101-102              V82 103-104               
   V83 105-105              V84 106-106              V85 107-108               
   V86 109-110              V87 111-112              V88 113-114               
   V89 115-115              V90 116-117              V91 118-119               
   V101 120-120             V102 121-122             V103 123-124              
   V104 125-125             V105 126-127             V106 128-129              
   V107 130-130             V108 131-131             V109 132-133              
   V110 134-135             V111 136-136             V112 137-138              
   V113 139-140             V114 141-141             V115 142-142              
   V116 143-144             V117 145-146             V118 147-147              
   V119 148-149             V120 150-151             V121 152-152              
   V122 153-154             V123 155-156             V124 157-157              
   V125 158-159             V126 160-161             V127 162-162              
   V128 163-164             V129 165-166             V130 167-168              
   V131 169-170             V132 171-171             V133 172-173              
   V134 174-175             V135 176-176             V136 177-177              
   V137 178-179             V138 180-181             V139 182-182              
   V140 183-183             V141 184-184             V142 185-185              
   V143 186-187             V144 188-189             V145 190-191              
   V146 192-193             V147 194-194             V148 195-195              
   V149 196-196             V150 197-197             V151 198-198              
   V152 199-199             V153 200-200             V154 201-201              
   V155 202-202             V156 203-203             V157 204-204              
   V158 205-205             V159 206-206             V160 207-207              
   V161 208-208             V162 209-209             V163 210-210              
   V164 211-211             V165 212-212             V166 213-213              
   V167 214-214             V168 215-215             V169 216-216              
   V170 217-217             V171 218-218             V172 219-219              
   V173 220-220             V174 221-221             V175 222-222              
   V176 223-223             V177 224-224             V178 225-225              
   V179 226-226             V180 227-227             V181 228-228              
   V182 229-229             V183 230-230             V184 231-232              
   V185 233-234             V186 235-235             V187 236-236              
   V188 237-238             V189 239-240             V190 241-241              
   V191 242-242             V192 243-243             V193 244-245              
   V194 246-247             V195 248-248             V196 249-250              
   V197 251-252             V198 253-253             V201 254-254              
   V202 255-255             V203 256-257             V204 258-258              
   V205 259-262             V206 263-266             V207 267-268              
   V208 269-269             V209 270-273             V210 274-277              
   V211 278-279             V212 280-280             V213 281-284              
   V214 285-288             V215 289-290             V216 291-291              
   V217 292-295             V218 296-299             V219 300-301              
   V220 302-302             V221 303-306             V222 307-310              
   V223 311-311             V224 312-312             V225 313-314              
   V226 315-315             V227 316-316             V228 317-317              
   V229 318-318             V230 319-320             V231 321-322              
   V232 323-323             V233 324-325             V234 326-327              
   V235 328-328             V236 329-332             V237 333-334              
   V238 335-336             V239 337-337             V240 338-339              
   V241 340-341             V242 342-342             V243 343-344              
   V244 345-346             V245 347-347             V246 348-351              
   V247 352-353             V248 354-355             V249 356-356              
   V250 357-358             V251 359-360             V252 361-361              
   V253 362-363             V254 364-365             V255 366-366              
   V256 367-367             V257 368-368             V258 369-369              
   V259 370-370             V260 371-371             V261 372-372              
   V262 373-373             V263 374-374             V264 375-375              
   V265 376-376             V266 377-377             V267 378-378              
   V268 379-379             V275 380-381             V276 382-383              
   V277 384-385             V278 386-387             V279 388-389              
   V280 390-391             V281 392-395             V282 396-396              
   V283 397-398             V284 399-400             V285 401-402              
   V286 403-404             V287 405-406             V288 407-408              
   V289 409-412             V301 413-414             V302 415-416              
   V303 417-417             V304 418-418             V305 419-422              
   V306 423-423             V307 424-424             V308 425-430              
   V309 431-433             V310 434-434             V311 435-440              
   V312 441-441             V313 442-442             V314 443-444              
   V315 445-445             V316 446-446             V317 447-452              
   V318 453-456             V319 457-457             V320 458-458              
   V321 459-460             V322 461-462             V323 463-464              
   V324 465-468             V325 469-469             V326 470-475              
   V327 476-479             V328 480-480             V329 481-481              
   V330 482-483             V331 484-485             V332 486-487              
   V333 488-491             V334 492-492             V335 493-493              
   V336 494-494             V337 495-495             V338 496-496              
   V339 497-498             V351 499-499             V352 500-501              
   V353 502-503             V354 504-508             V355 509-509              
   V356 510-510             V357 511-511             V358 512-516              
   V359 517-517             V360 518-518             V361 519-522              
   V362 523-523             V363 524-525             V364 526-526              
   V365 527-527             V366 528-528             V367 529-530              
   V368 531-532             V369 533-537             V370 538-538              
   V371 539-539             V372 540-540             V373 541-545              
   V374 546-546             V375 547-547             V376 548-551              
   V377 552-552             V378 553-554             V379 555-555              
   V380 556-556             V381 557-557             V382 558-559              
   V383 560-561             V384 562-566             V385 567-567              
   V386 568-568             V387 569-569             V388 570-574              
   V389 575-575             V390 576-576             V391 577-580              
   V392 581-581             V393 582-583             V394 584-584              
   V395 585-585             V425 586-586             V426 587-587              
   V427 588-589             V428 590-591             V429 592-593              
   V430 594-594             V431 595-597             V432 598-602              
   V433 603-603             V434 604-604             V435 605-605              
   V436 606-606             V437 607-607             V438 608-612              
   V439 613-616             V440 617-617             V441 618-618              
   V442 619-620             V443 621-621             V444 622-625              
   V445 626-626             V446 627-627             V447 628-628              
   V448 629-630             V449 631-632             V450 633-634              
   V451 635-635             V452 636-638             V453 639-643              
   V454 644-644             V455 645-645             V456 646-646              
   V457 647-647             V458 648-648             V459 649-653              
   V460 654-657             V461 658-658             V462 659-659              
   V463 660-661             V464 662-662             V465 663-666              
   V466 667-667             V467 668-668             V468 669-669              
   V469 670-671             V470 672-673             V471 674-675              
   V472 676-676             V473 677-679             V474 680-684              
   V475 685-685             V476 686-686             V477 687-687              
   V478 688-688             V479 689-689             V480 690-694              
   V481 695-698             V482 699-699             V483 700-700              
   V484 701-702             V485 703-703             V486 704-707              
   V487 708-708             V488 709-709             V501 710-710              
   V502 711-712             V503 713-714             V504 715-716              
   V505 717-721             V506 722-722             V507 723-723              
   V508 724-724             V509 725-729             V510 730-730              
   V511 731-734             V512 735-735             V513 736-736              
   V514 737-738             V515 739-739             V516 740-740              
   V517 741-741             V518 742-743             V519 744-745              
   V520 746-747             V521 748-752             V522 753-753              
   V523 754-754             V524 755-755             V525 756-760              
   V526 761-761             V527 762-765             V528 766-766              
   V529 767-767             V530 768-769             V531 770-770              
   V532 771-771             V533 772-772             V534 773-774              
   V535 775-776             V536 777-778             V537 779-783              
   V538 784-784             V539 785-785             V540 786-786              
   V541 787-791             V542 792-792             V543 793-796              
   V544 797-797             V545 798-798             V546 799-800              
   V547 801-801             V548 802-802             V601 803-803              
   V602 804-805             V603 806-807             V604 808-809              
   V605 810-814             V606 815-815             V607 816-816              
   V608 817-817             V609 818-822             V610 823-823              
   V611 824-827             V612 828-828             V613 829-829              
   V614 830-831             V615 832-832             V616 833-833              
   V617 834-834             V618 835-836             V619 837-838              
   V620 839-840             V621 841-845             V622 846-846              
   V623 847-847             V624 848-848             V625 849-853              
   V626 854-854             V627 855-858             V628 859-859              
   V629 860-860             V630 861-862             V631 863-863              
   V632 864-864             V633 865-865             V634 866-867              
   V635 868-869             V636 870-871             V637 872-876              
   V638 877-877             V639 878-878             V640 879-879              
   V641 880-884             V642 885-885             V643 886-889              
   V644 890-890             V645 891-891             V646 892-893              
   V647 894-894             V648 895-895             V701 896-896              
   V702 897-898             V703 899-902             V704 903-903              
   V705 904-908             V706 909-910             V707 911-912              
   V708 913-914             V709 915-915             V710 916-916              
   V711 917-917             V712 918-918             V713 919-919              
   V714 920-921             V715 922-925             V716 926-926              
   V717 927-931             V718 932-933             V719 934-935              
   V720 936-937             V721 938-938             V722 939-939              
   V723 940-940             V724 941-941             V725 942-942              
   V726 943-944             V727 945-948             V728 949-949              
   V729 950-954             V730 955-956             V731 957-958              
   V732 959-960             V733 961-961             V734 962-962              
   V735 963-963             V736 964-964             V737 965-966              
   V738 967-967             V739 968-972             V740 973-974              
   V741 975-975             V742 976-977             V743 978-978              
   V744 979-983             V745 984-985             V746 986-986              
   V747 987-988             V748 989-989             V749 990-994              
   V750 995-996             V751 997-997             V752 998-998              
   V753 999-1000            V754 1001-1002           V801 1003-1003            
   V802 1004-1005           V803 1006-1007           V804 1008-1008            
   V805 1009-1009           V806 1010-1010           V807 1011-1014            
   V808 1015-1015           V809 1016-1016           V810 1017-1017            
   V811 1018-1019           V812 1020-1021           V813 1022-1023            
   V814 1024-1025           V815 1026-1027           V816 1028-1029            
   V817 1030-1030           V818 1031-1032           V819 1033-1033            
   V820 1034-1034           V821 1035-1036           V822 1037-1038            
   V823 1039-1040           V824 1041-1042           V825 1043-1043            
   V826 1044-1044           V827 1045-1045           V828 1046-1046            
   V829 1047-1047           V830 1048-1048           V831 1049-1050            
   V832 1051-1052           V833 1053-1053           V834 1054-1054            
   V835 1055-1055           V836 1056-1056           V837 1057-1057            
   V838 1058-1058           V839 1059-1059           V840 1060-1060            
   V841 1061-1061           V842 1062-1062           V843 1063-1063            
   V844 1064-1064           V845 1065-1065           V846 1066-1067            
   V847 1068-1069           V848 1070-1070           V849 1071-1072            
   V850 1073-1074           V851 1075-1075           V852 1076-1076            
   V853 1077-1077           V854 1078-1079           V855 1080-1080            
   V856 1081-1081           V857 1082-1082           V858 1083-1083            
   V859 1084-1084           V860 1085-1085           V861 1086-1087            
   V862 1088-1089           V863 1090-1090           V864 1091-1091            
   V865 1092-1093           V866 1094-1095           V867 1096-1097            
   V868 1098-1099           V869 1100-1101           V870 1102-1103            
   V871 1104-1104           V872 1105-1105           V873 1106-1106            
   V874 1107-1107           V875 1108-1108           V876 1109-1110            
   V877 1111-1111           V901 1112-1112           V902 1113-1113            
   V903 1114-1114           V904 1115-1115           V905 1116-1117            
   V906 1118-1119           V907 1120-1120           V908 1121-1122            
   V909 1123-1124           V910 1125-1125           V911 1126-1126            
   V912 1127-1127           V913 1128-1128           V914 1129-1130            
   V915 1131-1132           V916 1133-1133           V917 1134-1134            
   V918 1135-1136           V919 1137-1137           V920 1138-1139            
   V921 1140-1140           V922 1141-1142           V923 1143-1143            
   V924 1144-1145           V925 1146-1146           V926 1147-1148            
   V927 1149-1149           V928 1150-1151           V929 1152-1152            
   V930 1153-1154           V931 1155-1155           V932 1156-1157            
   V933 1158-1158           V934 1159-1160           V935 1161-1165            
   V936 1166-1170           V937 1171-1172           V938 1173-1177            
   V939 1178-1182           V940 1183-1184           V941 1185-1189            
   V942 1190-1194           V943 1195-1195           V944 1196-1196            
   V945 1197-1197           V946 1198-1199           V947 1200-1201            
   V948 1202-1203           V949 1204-1204           V950 1205-1205            
   V951 1206-1206           V952 1207-1208           V953 1209-1210            
   V954 1211-1212           V955 1213-1214           V956 1215-1216            
   V957 1217-1218           V958 1219-1220           V959 1221-1221            
   V960 1222-1222           V961 1223-1223           V962 1224-1224            
   V963 1225-1226           V964 1227-1227           V965 1228-1228            
   V966 1229-1229           V967 1230-1230           V968 1231-1231            
   V969 1232-1233           V970 1234-1234           V1001 1235-1236           
   V1002 1237-1238          V1003 1239-1240          V1004 1241-1245           
   V1005 1246-1246          V1006 1247-1247          V1007 1248-1248           
   V1008 1249-1253          V1009 1254-1254          V1010 1255-1258           
   V1011 1259-1259          V1012 1260-1260          V1013 1261-1262           
   V1014 1263-1263          V1015 1264-1264          V1021 1265-1265           
   V1022 1266-1266          V1023 1267-1267          V1024 1268-1268           
   V1025 1269-1269          V1026 1270-1270          V1027 1271-1272           
   V1028 1273-1273          V1029 1274-1275          V1030 1276-1276           
   V1031 1277-1282 .2       V1032 1283-1288 .2       V1033 1289-1290           
   V1034 1291-1292          V1035 1293-1293          V1036 1294-1294;          
                                                                               
* SAS LABEL STATEMENT;                                                         
                                                                               
LABEL                                                                          
   V1 = "RESPONDENT ID #"                                                      
   V3 = "A1:G/B BUY NSTLMT PLAN"                                               
   V4 = "A1A1:WHY G/B BUY NSTLMT"                                              
   V5 = "A1A2:WHY G/B BUY NSTLMT"                                              
   V6 = "A2A:ALRITE BRRW-VACATION"                                             
   V7 = "A2B:ALRITE BRRW-LVNG XPN"                                             
   V8 = "A2C:ALRITE BRRW-CNSL BLS"                                             
   V9 = "A2D:ALRITE BRRW-FUR,JWLS"                                             
   V10 = "A2E:ALRITE BRRW-BT,SNMBL"                                            
   V11 = "A2F:ALRITE BRRW-CAR"                                                 
   V12 = "A2G:ALRITE BRRW-ILLNESS"                                             
   V13 = "A2H:ALRITE BRRW-EDCTN XP"                                            
   V14 = "A2I:ALRITE BRRW-FURNITUR"                                            
   V15 = "A3A:LENDRS-FAIR TO CUST"                                             
   V16 = "A3B:LENDRS-GD IF NT ARND"                                            
   V17 = "A3C:LENDRS-MST PPL SAT"                                              
   V18 = "A3D:LENDRS-PRVD USFL SVC"                                            
   V19 = "A3E:LENDRS-TAKE ADVNTAGE"                                            
   V20 = "A4:DIFFICLT SAV ADVNC $"                                             
   V21 = "A5-1:WHT WNT KNW CRD SHP"                                            
   V22 = "A5-2:WHT WNT KNW CRD SHP"                                            
   V23 = "A5-3:WHT WNT KNW CRD SHP"                                            
   V24 = "A5-4:WHT WNT KNW CRD SHP"                                            
   V25 = "A5A:MST MPTNT KNW CRD SH"                                            
   V26 = "A6A:CAR-RANK    AMT FNCD"                                            
   V27 = "A6A-A:DRBL-RANK-AMT FNCD"                                            
   V28 = "A6B-A:H AND R-RANK- AMT FNCD"                                        
   V29 = "A6B:CAR-RANK-   FNC CHG"                                             
   V30 = "A6A-B:DRBL-RANK-FNC CHG"                                             
   V31 = "A6B-B:H AND R-RANK- FNC CHG"                                         
   V32 = "A6-C:CAR-RANK-  MTHL PMT"                                            
   V33 = "A6A-C:DRBL-RANK-MTHL PMT"                                            
   V34 = "A6B-C:H AND R-RANK- MTHL PMT"                                        
   V35 = "A6-D:CAR-RANK-  %RT INTR"                                            
   V36 = "A6A-D:DRBL-RANK-%RT INTR"                                            
   V37 = "A6B-D:H AND R-RANK- %RT INTR"                                        
   V38 = "A6E:CAR-RANK-   LATE CHG"                                            
   V39 = "A6A-E:DRBL-RANK-LATE CHG"                                            
   V40 = "A6B-E:H AND R-RANK- LATE CHG"                                        
   V41 = "A6-F:CAR-RANK-  COLLATRL"                                            
   V42 = "A6A-F:DRBL-RANK-COLLATRL"                                            
   V43 = "A6B-F:H AND R-RANK- COLLATRL"                                        
   V44 = "A6-G:CAR-RANK-  ERL PYOF"                                            
   V45 = "A6A-G:DRBL-RANK-ERL PYOF"                                            
   V46 = "A6B-G:H AND R-RANK- ERL PYOF"                                        
   V47 = "A7-#1:RCMD FRND-STOR/DLR"                                            
   V48 = "A7-#2:RCMD FRND-BANK"                                                
   V49 = "A7-#3:RCMD FRND-FINC CO"                                             
   V50 = "A7-#4:RCMD FRND-CRDT UNI"                                            
   V51 = "A7A(1):WHY REC LOAN SRCE"                                            
   V52 = "A7A(2):WHY REC LOAN SRCE"                                            
   V53 = "A7A(3):WHY REC LOAN SRCE"                                            
   V54 = "A8:CAUTN FRND-INST N LST"                                            
   V55 = "A8A-#1:CAUTN FRN-STR/DLR"                                            
   V56 = "A8A-#2:CAUTN FRN-BANK"                                               
   V57 = "A8A-#3:CAUTN FRN-FINC CO"                                            
   V58 = "A8A-#4:CAUTN FRN-CRDT UN"                                            
   V59 = "A8B(1):WHY CAUTN FRND"                                               
   V60 = "A8B(2):WHY CAUTN FRND"                                               
   V61 = "A9A:CRDT DFCLT-STR/DLR"                                              
   V62 = "A9B:CRDT DFCLT-BANK"                                                 
   V63 = "A9C:CRDT DFCLT-FINC CPY"                                             
   V64 = "A9D:CRDT DFCLT-CRDT UNI"                                             
   V65 = "A10:CRDT DFCLT GET SELF"                                             
   V66 = "A10A-#1:CRD DFC SLF-STR"                                             
   V67 = "A10A-#2:CRD DFC SLF-BANK"                                            
   V68 = "A10A-#3:CRD DFC SLF-FINC"                                            
   V69 = "A10A-#4:CRD DFC SLF-CR.U"                                            
   V70 = "A10B(1):WHY DFCLT SELF"                                              
   V71 = "A10B(2):WHY DFCLT SELF"                                              
   V72 = "A11:COST $1000 LOAN"                                                 
   V73 = "A11A:PROBE USED IN A11?"                                             
   V74 = "A12:YRLY NTRST ON $1000"                                             
   V75 = "A12A:PROBE USED IN A12?"                                             
   V76 = "A13A:CRDT COST-STOR/DLR"                                             
   V77 = "A13B:CRDT COST-BANK"                                                 
   V78 = "A13C:CRDT COST-FINC CPNY"                                            
   V79 = "A13D:CRDT COST-CRDT UNI"                                             
   V80 = "A14:CREDITORS GV GD INFO"                                            
   V81 = "A14A:WHAT LOAN INF HLPFL"                                            
   V82 = "A14B:WHAT LOAN INF HLPFL"                                            
   V83 = "A15:CRDT INF HARD TO GET"                                            
   V84 = "A16:LOAN THRU DEALER G/B"                                            
   V85 = "A17(1):DEALR LOAN-ADVANT"                                            
   V86 = "A17(2):DEALR LOAN-ADVANT"                                            
   V87 = "A18(1):DEALR LOAN-DISADV"                                            
   V88 = "A18(2):DEALR LOAN-DISADV"                                            
   V89 = "A19:CNSMR PRTCTN LWS EFF"                                            
   V90 = "A19A1:HOW CHNG CNSM LAWS"                                            
   V91 = "A19A2:HOW CHNG CNSM LAWS"                                            
   V101 = "B1:R TRTD NFRLY CRD TRNS"                                           
   V102 = "B1A-#1A:WHAT WAS CRD PBM"                                           
   V103 = "B1A-#1B:WHAT WAS CRD PBM"                                           
   V104 = "B1B-#1:TRY FIX CRDT PBM"                                            
   V105 = "B1C-#1A:WHT DO FX CRD PB"                                           
   V106 = "B1C-#1B:WHT DO FX CRD PB"                                           
   V107 = "B1D-#1:PBM CORR TO STFCN"                                           
   V108 = "B1E:OTHR UNFR CRDT PBMS"                                            
   V109 = "B1A-#2A:WHT OTHR CRD PBM"                                           
   V110 = "B1A-#2B:WHT OTHR CRD PBM"                                           
   V111 = "B1B-#2:TRY FIX OTR CR PB"                                           
   V112 = "B1C-#2A:WHT DO FX OTR PB"                                           
   V113 = "B1C-#2B:WHT DO FX OTR PB"                                           
   V114 = "B1D-#2:OTR PBM CORR STFC"                                           
   V115 = "B2:NFR CRD PRC WNT CHNGD"                                           
   V116 = "B2A-#1A:WHT PRCT CHANGE"                                            
   V117 = "B2A-#1B:WHT PRCT CHANGE"                                            
   V118 = "B2B-#1:ALL CRDTRS SAME"                                             
   V119 = "B2C-#1A:DIFFRNCS-CRDTRS"                                            
   V120 = "B2C-#1B:DIFFRNCS-CRDTRS"                                            
   V121 = "B2D:OTR PRTCS WNT CHNGD"                                            
   V122 = "B2A-#2A:WHT OTR PRCT CHD"                                           
   V123 = "B2A-#2B:WHT OTR PRCT CHD"                                           
   V124 = "B2B-#2:ALL CRDTRS SM/DIF"                                           
   V125 = "B2C-#2A:HOW CRDTRS DIFF"                                            
   V126 = "B2C-#2B:HOW CRDTRS DIFF"                                            
   V127 = "B3:COMPLAIN PST CRDT EXP"                                           
   V128 = "B3A1:WHO CMPLN PST EXPS"                                            
   V129 = "B3A2:WHO CMPLN PST EXPS"                                            
   V130 = "B3B1:WHT CAUSD CMPLNT"                                              
   V131 = "B3B2:WHT CAUSD CMPLNT"                                              
   V132 = "B3C:OUTCOME OF COMPLAINT"                                           
   V133 = "B4(1):AGNCY CMPLN CRD PB"                                           
   V134 = "B4(2):AGNCY CMPLN CRD PB"                                           
   V135 = "B4A:PROBE USED?"                                                    
   V136 = "B5:OTR AGNCY CMPL CRD PB"                                           
   V137 = "B5A1:OTR AGNC CMP CRD PB"                                           
   V138 = "B5A2:OTR AGNC CMP CRD PB"                                           
   V139 = "B6A:READ TRTH.LNDNG STMT"                                           
   V140 = "B6B:TRTH.LNDNG STMT COMP"                                           
   V141 = "B6C:TRTH.LNDNG.S NT USFL"                                           
   V142 = "B6D:TRTH.LNDNG=+CDT CONF"                                           
   V143 = "B7(1):INF CDTRS US MK LN"                                           
   V144 = "B7(2):INF CDTRS US MK LN"                                           
   V145 = "B7(3):INF CDTRS US MK LN"                                           
   V146 = "B7(4):INF CDTRS US MK LN"                                           
   V147 = "B8A:CRDT FACT-TIME N JOB"                                           
   V148 = "B9A:WNT ILLGL-TIME N JOB"                                           
   V149 = "B10A:ARE ILGL-TIME N JOB"                                           
   V150 = "B8B:CRDT FACT-LN PR ADRS"                                           
   V151 = "B9B:WNT ILLGL-LN PR ADRS"                                           
   V152 = "B10B:ARE ILGL-LN PR ADRS"                                           
   V153 = "B8C:CRDT FACT-RACE"                                                 
   V154 = "B9C:WNT ILLGL-RACE"                                                 
   V155 = "B10C:ARE ILGL-RACE"                                                 
   V156 = "B8D:CRDT FACT-HV CHK ACT"                                           
   V157 = "B9D:WNT ILLGL-HV CHK ACT"                                           
   V158 = "B10D:ARE ILGL-HV CHK ACT"                                           
   V159 = "B8E:CRDT FACT-OWN HOME"                                             
   V160 = "B9E:WNT ILLGL-OWN HOME"                                             
   V161 = "B10E:ARE ILGL-OWN HOME"                                             
   V162 = "B8F:CRDT FACT-OTR M.PMTS"                                           
   V163 = "B9F:WNT ILLGL-OTR M.PMTS"                                           
   V164 = "B10F:ARE ILGL-OTR M.PMTS"                                           
   V165 = "B8G:CRDT FACT-AGE"                                                  
   V166 = "B9G:WNT ILLGL-AGE"                                                  
   V167 = "B10G:ARE ILGL-AGE"                                                  
   V168 = "B8H:CRDT FACT-SEX"                                                  
   V169 = "B9H:WNT ILLGL-SEX"                                                  
   V170 = "B10H:ARE ILGL-SEX"                                                  
   V171 = "B8I:CRDT FACT-INCOME"                                               
   V172 = "B9I:WNT ILLGL-INCOME"                                               
   V173 = "B10I:ARE ILGL-INCOME"                                               
   V174 = "B8J:CRDT FACT-MARITAL ST"                                           
   V175 = "B9J:WNT ILLGL-MARITAL ST"                                           
   V176 = "B10J:ARE ILGL-MARITAL ST"                                           
   V177 = "B8K:CRDT FACT-FMLY SIZE"                                            
   V178 = "B9K:WNT ILLGL-FMLY SIZE"                                            
   V179 = "B10K:ARE ILGL-FMLY SIZE"                                            
   V180 = "B8L:CRDT FACT-PRV CRD EX"                                           
   V181 = "B9L:WNT ILLGL-PRV CRD EX"                                           
   V182 = "B10L:ARE ILGL-PRV CRD EX"                                           
   V183 = "B11:TRN DWN F.CRD LST YR"                                           
   V184 = "B11A1:REASON NO CREDIT"                                             
   V185 = "B11A2:REASON NO CREDIT"                                             
   V186 = "B11B:WAS RSN FR N.CRD OK"                                           
   V187 = "B12:UNABL GT ENUF CRDT"                                             
   V188 = "B12A1:WHY NOT ENUF CRDT"                                            
   V189 = "B12A2:WHY NOT ENUF CRDT"                                            
   V190 = "B12B:EXPL NT ENUF CRD OK"                                           
   V191 = "B13:CHECKPOINT"                                                     
   V192 = "B14:TRN DWN=AGE,SEX,RACE"                                           
   V193 = "B14A1:WHAT HAPPENED"                                                
   V194 = "B14A2:WHAT HAPPENED"                                                
   V195 = "B14B:TRY DO ANYTHING"                                               
   V196 = "B14C1:WHAT TRY DO-TRN DN"                                           
   V197 = "B14C2:WHAT TRY DO-TRN DN"                                           
   V198 = "B14D:PBM CORRECTION OK"                                             
   V201 = "C1:G/B CREDIT CARDS"                                                
   V202 = "C2:R/FMLY HAVE CRDT CRDS"                                           
   V203 = "C3A:NUMBER GASOLINE CRDS"                                           
   V204 = "C4A:USE GASOLINE CARDS"                                             
   V205 = "C5A:AMT CHGD LST MO-GAS"                                            
   V206 = "C6A:BALANCE-GAS CARDS"                                              
   V207 = "C3B:NUMBR BANK CRDT CRDS"                                           
   V208 = "C4B:USE BANK CARDS"                                                 
   V209 = "C5B:AMT CHGD LST MO-BANK"                                           
   V210 = "C6B:BALANCE-BANK CARDS"                                             
   V211 = "C3C:NMBR GENL PURPS CRDS"                                           
   V212 = "C4C:USE GEN PURPOSE CRDS"                                           
   V213 = "C5C:AMT CHGD LST MO-GP"                                             
   V214 = "C6C:BALANCE-GEN PURP CRD"                                           
   V215 = "C3D:NMBR STORE CRDS/ACCT"                                           
   V216 = "C4D:USE STORE CARD/ACCT"                                            
   V217 = "C5D:AMT CHGD LST MO-STOR"                                           
   V218 = "C6D:BALANCE-STORE CRDS"                                             
   V219 = "C3E:NMBR OTHER CRDT CRDS"                                           
   V220 = "C4E:USE OTHER CRDT CARDS"                                           
   V221 = "C5E:AMT CHGD LST MO-OTHR"                                           
   V222 = "C6E:BALANCE-OTHR CRDS"                                              
   V223 = "C7:CHECKPOINT-HAVE CARDS"                                           
   V224 = "C8:ALWAYS PAY IN FULL"                                              
   V225 = "C9:%TAGE RATE ON BALANCE"                                           
   V226 = "C10:<3%)-MONTHLY/ANNUAL"                                            
   V227 = "C11:CHECK STMT W/SLIPS"                                             
   V228 = "C12:GET DSCRPTN/SLIPS"                                              
   V229 = "C12A:DSCRPTNS ADEQUATE"                                             
   V230 = "C12B1:HW CHNG DSCR MK AD"                                           
   V231 = "C12B2:HW CHNG DSCR MK AD"                                           
   V232 = "C13:CARD BILL ERRS LST.Y"                                           
   V233 = "C13A1-1:WHT BILL ERR MD"                                            
   V234 = "C13A1-2:WHT BILL ERR MD"                                            
   V235 = "C13B1:BILL ERR-WHICH CRD"                                           
   V236 = "C13C1:SIZE BILLNG ERR-$"                                            
   V237 = "C13D1-1:WHT DO ABT ERROR"                                           
   V238 = "C13D1-2:WHT DO ABT ERROR"                                           
   V239 = "C13E1:BLLNG ERR CORR OK"                                            
   V240 = "C13F1-1:PBLM IN FXNG ERR"                                           
   V241 = "C13F1-2:PBLM IN FXNG ERR"                                           
   V242 = "C13G:OTHR ERR-CRD CD BLS"                                           
   V243 = "C13A2-1:WHT OTH BIL ERRS"                                           
   V244 = "C13A2-2:WHT OTH BIL ERRS"                                           
   V245 = "C13B2:OTR BIL ERR-WCH CD"                                           
   V246 = "C13C2:SIZE OTR BLNG ERR$"                                           
   V247 = "C13D2-1:WHT DO OTHR ERR"                                            
   V248 = "C13D2-2:WHT DO OTHR ERR"                                            
   V249 = "C13E2:OTR BLNG ERR CORR"                                            
   V250 = "C13F2-1:PBM IN FX OTR ER"                                           
   V251 = "C13F2-2:PBM IN FX OTR ER"                                           
   V252 = "C14:BILLNG PRCTS CHANGE"                                            
   V253 = "C14A1:WHT BLNG PRCTS CHN"                                           
   V254 = "C14A2:WHT BLNG PRCTS CHN"                                           
   V255 = "C15:FDRL LGSL-CRD BL ERR"                                           
   V256 = "C16:MARITAL STATUS OF R"                                            
   V257 = "C16A:R MARRIED BEFORE"                                              
   V258 = "C17:CHECKPOINT-SEX OF R"                                            
   V259 = "C18:CRD CRDS BFR MARRIAG"                                           
   V260 = "C18A1:PRMAR C.CD-CANCEL"                                            
   V261 = "C18A2:PRMAR C.CD-CHN NAM"                                           
   V262 = "C18A3:PRMAR C.CD-REAPPLY"                                           
   V263 = "C18B:CHNG NAME/NEW ACCNT"                                           
   V264 = "C19:CRD CRDS BFR MARRIAG"                                           
   V265 = "C19A1:PRMAR C.CD-CANCEL"                                            
   V266 = "C19A2:PRMAR C.CD-CHN NAM"                                           
   V267 = "C19A3:PRMAR C.CD-REAPPLY"                                           
   V268 = "C19B:CHNG NAME/NEW ACCNT"                                           
   V275 = "D1:EMPLOYMENT STATUS-R"                                             
   V276 = "D2:#HOURS WORKED/WEEK-R"                                            
   V277 = "D2A-D4C:OCCUPATN CODE -R"                                           
   V278 = "D7:# YEARS EDUCATION-R"                                             
   V279 = "D7-D7C:EDUCTN SUMMARY-R"                                            
   V280 = "D8:R BIRTHDATE-MONTH"                                               
   V281 = "D8:R BIRTHDATE-YEAR"                                                
   V282 = "D9:CHKPT-R HAVE SPOUSE"                                             
   V283 = "D10:EMPLOYMENT STATUS-SP"                                           
   V284 = "D11:#HOURS WORKED/WK-SP"                                            
   V285 = "D11A-D13C:OCC CODE-SPOU"                                            
   V286 = "D16:#YEARS EDUCATION-SP"                                            
   V287 = "D16-D16C:EDUCTN SUM-SPOU"                                           
   V288 = "D17:SPOUSE BIRTHDATE-MON"                                           
   V289 = "D17:SPOUSE BIRTHDATE-YR"                                            
   V301 = "F1:YEAR MOVED IN HOUSE"                                             
   V302 = "F2:LENGTH COUNTY RESIDNC"                                           
   V303 = "F3:HOUSING-OWN/RENT/OTHR"                                           
   V304 = "F4:HOUSING=OTHER)WHAT?"                                             
   V305 = "F5:RENT PER MONTH"                                                  
   V306 = "F5A:RENT INCLUDE UTILITY"                                           
   V307 = "F5B:RENT (UN)FURNISHED"                                             
   V308 = "F6:PRESENT VALUE-HOUSE"                                             
   V309 = "F6A:MOBL HM) SITE RENT"                                             
   V310 = "F7:MOBILE HM) NEW/USED"                                             
   V311 = "F8:ORIGINAL COST OF HOME"                                           
   V312 = "F9:MORTGAGE ON HOME"                                                
   V313 = "F10:2ND MORTGAGE ON HOME"                                           
   V314 = "F10A:YEAR 2ND MORTGAGE"                                             
   V315 = "F10B:WHAT 2ND MORTG FOR"                                            
   V316 = "F10B:WHAT 2ND MORTG FOR"                                            
   V317 = "F11(1):SIZE 1ST MORTGAGE"                                           
   V318 = "F12(1):$/MO 1ST MORTGAGE"                                           
   V319 = "F12A(1):PMT INCL TAX/INS"                                           
   V320 = "F13(1):PMT SAME NOW-THEN"                                           
   V321 = "F14(1):HOW PMTS CHANGED"                                            
   V322 = "F14(1):HOW PMTS CHANGED"                                            
   V323 = "F15(1):#YRS LEFT MORTGAG"                                           
   V324 = "F16(1):ANNUAL %RATE INTR"                                           
   V325 = "F17(1):PAY BANK/S AND L/OTHR"                                       
   V326 = "F11(2):SIZE 2ND MORTGAGE"                                           
   V327 = "F12(2):$/MO 2ND MORTGAGE"                                           
   V328 = "F12A(2):PMT INCL TAX/INS"                                           
   V329 = "F13(2):PMT SAME NOW-THEN"                                           
   V330 = "F14(2):HOW PMTS CHANGED"                                            
   V331 = "F14(2):HOW PMTS CHANGED"                                            
   V332 = "F15(2):#YRS LEFT MORTGAG"                                           
   V333 = "F16(2):ANNUAL %RATE INTR"                                           
   V334 = "F17(2):PAY BANK/S AND L/OTHR"                                       
   V335 = "F18:REFINANCE 1ST MORTG"                                            
   V336 = "F18A:PAY OFF/ADD TO OLD"                                            
   V337 = "F18B:WHY REFINANCE 1ST M"                                           
   V338 = "F18B:WHY REFINANCE 1ST M"                                           
   V339 = "F18C:WHEN REFINANCE-YEAR"                                           
   V351 = "G1:ADDITIONS AND REPAIRS7677"                                       
   V352 = "G3M-1:A AND R START MONTH"                                          
   V353 = "G3Y-1:A AND R START YEAR"                                           
   V354 = "G4 -1:$ COST"                                                       
   V355 = "G5 -1:CASH/CREDIT/CHARGE"                                           
   V356 = "G5A-1: CASH-LOAN/SAV/OTR"                                           
   V357 = "G6 -1:CHARGE-BANK/STORE"                                            
   V358 = "G7 -1:AMOUNT BORROWED"                                              
   V359 = "G8 -1:ANYTHING LEFT PAY"                                            
   V360 = "G9 -1:AMOUNT INCL MORTG"                                            
   V361 = "G10-1:AMT LOAN PAYMENTS"                                            
   V362 = "G10-1:PAYMENT FREQUENCY"                                            
   V363 = "G11-1:# PAYMENTS"                                                   
   V364 = "G12-1:PAY BNK/FCO/CU/STR"                                           
   V365 = "G13-1:STOR/CNTRCT MK ARR"                                           
   V366 = "G14-1:OTR A AND R 76/77"                                            
   V367 = "G3M-2:A AND R START MONTH"                                          
   V368 = "G3Y-2:A AND R START YEAR"                                           
   V369 = "G4 -2:$COST"                                                        
   V370 = "G5 -2:CASH/CREDIT/CHARGE"                                           
   V371 = "G5A-2:CASH-LOAN/SAV/OTHR"                                           
   V372 = "G6 -2:CHARGE-BANK/STORE"                                            
   V373 = "G7 -2:AMOUNT BORROWED"                                              
   V374 = "G8 -2:ANYTHING LEFT PAY"                                            
   V375 = "G9 -2:AMOUNT INCL MORTG"                                            
   V376 = "G10-2:AMT LOAN PAYMENTS"                                            
   V377 = "G10-2:PAYMENT FREQUENCY"                                            
   V378 = "G11-2:# PAYMENTS"                                                   
   V379 = "G12-2:PAY BNK/FCO/CU/STR"                                           
   V380 = "G13-2:STOR/CNTRCT MK ARR"                                           
   V381 = "G14-2:OTR A AND R 76/77"                                            
   V382 = "G3M-3:A AND R START MONTH"                                          
   V383 = "G3Y-3:A AND R START YEAR"                                           
   V384 = "G4 -3:$COST"                                                        
   V385 = "G5 -3:CASH/CREDIT/CHARGE"                                           
   V386 = "G5A-3:CASH-LOAN/SAV/OTHR"                                           
   V387 = "G6 -3:CHARGE-BANK/STORE"                                            
   V388 = "G7 -3:AMOUNT BORROWED"                                              
   V389 = "G8 -3:ANYTHING LEFT PAY"                                            
   V390 = "G9 -3:AMOUNT INCL MORTG"                                            
   V391 = "G10-3:AMT LOAN PAYMENTS"                                            
   V392 = "G10-3:PAYMENT FREQUENCY"                                            
   V393 = "G11-3:# PAYMENTS"                                                   
   V394 = "G12-3:PAY BNK/FCO/CU/STR"                                           
   V395 = "G13-3:STOR/CNTRCT MK ARR"                                           
   V425 = "H1:FU OWN/LEASE CAR"                                                
   V426 = "H2:# VEHICLES FU OWN/LS"                                            
   V427 = "H3-1:YEAR BUY VEHICLE #1"                                           
   V428 = "H4M1:WHN BUY/LEASE CAR-M"                                           
   V429 = "H4Y1:WHN BUY/LEASE CAR-Y"                                           
   V430 = "H4A1:BOUGHT NEW/USED-V#1"                                           
   V431 = "H5-1:MAKE AND MODEL-VEH#1"                                          
   V432 = "H6-1:COST VEHICLE #1"                                               
   V433 = "H7-1:BUY CASH/CRDT/LEASD"                                           
   V434 = "H7A1:CASH-LOAN/SAV/OTHR"                                            
   V435 = "H7B1:MAKING PAYMENTS NOW"                                           
   V436 = "H8-1:MAKING PAYMENTS NOW"                                           
   V437 = "H9-1:LEASING CONSIDERED"                                            
   V438 = "H10-1:SIZE OF LOAN-VEH#1"                                           
   V439 = "H11-1:AMT LOAN PAYMENTS"                                            
   V440 = "H11-1:PAYMENT FREQUENCY"                                            
   V441 = "H12-1:PMT NCL CAR/CR INS"                                           
   V442 = "H13-1:#PAYMENTS LOAN V#1"                                           
   V443 = "H14-1:FINAL PAYMENT SAME"                                           
   V444 = "H15-1:AMT FINAL PAYMENT"                                            
   V445 = "H16-1:PAY BNK/FCO/CU/STR"                                           
   V446 = "H17-1:DEALER MK ARR LOAN"                                           
   V447 = "H18-1:CHECKPOINT: >1 CAR"                                           
   V448 = "H3-2:YEAR BUY VEHICLE #2"                                           
   V449 = "H4M2:WHN BUY/LEASE CAR-M"                                           
   V450 = "H4Y2:WHN BUY/LEASE CAR-Y"                                           
   V451 = "H4A2:BOUGHT NEW/USED-V#2"                                           
   V452 = "H5-2:MAKE AND MODEL-VEH#2"                                          
   V453 = "H6-2:COST VEHICLE #2"                                               
   V454 = "H7-2:BUY CASH/CRDT/LEASD"                                           
   V455 = "H7A2:CASH-LOAN/SAV/OTHR"                                            
   V456 = "H7B2:MAKING PAYMENTS NOW"                                           
   V457 = "H8-2:MAKING PAYMENTS NOW"                                           
   V458 = "H9-2:LEASING CONSIDERED"                                            
   V459 = "H10-2:SIZE OF LOAN-VEH#2"                                           
   V460 = "H11-2:AMT LOAN PAYMENTS"                                            
   V461 = "H11-2:PAYMENT FREQUENCY"                                            
   V462 = "H12-2:PMT NCL CAR/CR INS"                                           
   V463 = "H13-2:#PAYMENTS LOAN V#2"                                           
   V464 = "H14-2:FINAL PAYMENT SAME"                                           
   V465 = "H15-2:AMT FINAL PAYMENT"                                            
   V466 = "H16-2:PAY BNK/FCO/CU/STR"                                           
   V467 = "H17-2:DEALER MK ARR LOAN"                                           
   V468 = "H19-2:CHECKPOINT: >2CARS"                                           
   V469 = "H3-3:YEAR BUY VEHICLE #3"                                           
   V470 = "H4M3:WHN BUY/LEASE CAR-M"                                           
   V471 = "H4Y3:WHN BUY/LEASE CAR-Y"                                           
   V472 = "H4A3:BOUGHT NEW/USED-V#3"                                           
   V473 = "H5-3:MAKE AND MODEL-VEH#3"                                          
   V474 = "H6-3:COST VEHICLE #3"                                               
   V475 = "H7-3:BUY CASH/CRDT/LEASD"                                           
   V476 = "H7A3:CASH-LOAN/SAV/OTHR"                                            
   V477 = "H7B3:MAKING PAYMENTS NOW"                                           
   V478 = "H8-3:MAKING PAYMENTS NOW"                                           
   V479 = "H9-3:LEASING CONSIDERED"                                            
   V480 = "H10-3:SIZE OF LOAN-VEH#3"                                           
   V481 = "H11-3:AMT LOAN PAYMENTS"                                            
   V482 = "H11-3:PAYMENT FREQUENCY"                                            
   V483 = "H12-3:PMT NCL CAR/CR INS"                                           
   V484 = "H13-3:#PAYMENTS LOAN V#3"                                           
   V485 = "H14-3:FINAL PAYMENT SAME"                                           
   V486 = "H15-3:AMT FINAL PAYMENT"                                            
   V487 = "H16-3:PAY BNK/FCO/CU/STR"                                           
   V488 = "H17-3:DEALER MK ARR LOAN"                                           
   V501 = "J1:BUY HH DURABLES-76/77"                                           
   V502 = "J2(1):HH DRBL-ITEM CODE"                                            
   V503 = "J3(1):HH DRBL-MNTH PURCH"                                           
   V504 = "J3(1):HH DRBL-YEAR PURCH"                                           
   V505 = "J4(1):COST OF ITEM"                                                 
   V506 = "J5(1):CASH/CRDT/CHG ACCT"                                           
   V507 = "J5A(1):CASH FRM LOAN/SAV"                                           
   V508 = "J6(1):BANKCARD/STOR ACCT"                                           
   V509 = "J7(1):AMOUNT FINANCED"                                              
   V510 = "J8(1):MAKING PAYMENT NOW"                                           
   V511 = "J9(1):PAYMENT AMOUNT"                                               
   V512 = "J9(1):PAYMENT FREQUENCY"                                            
   V513 = "J10(1):CRDT INSURN W/PMT"                                           
   V514 = "J11(1):#PMTS AGREED ON"                                             
   V515 = "J12(1):PAY BNK/FC/CU/STR"                                           
   V516 = "J13(1):CREDIT THRU STORE"                                           
   V517 = "J14(1):BUY OTHR HHDRBLS"                                            
   V518 = "J2(2):ITEM CODE"                                                    
   V519 = "J3(2):WHEN BOUGHT-MONTH"                                            
   V520 = "J3(2):WHEN BOUGHT-YEAR"                                             
   V521 = "J4(2):COST OF ITEM"                                                 
   V522 = "J5(2):CASH/CRDT/CHG ACCT"                                           
   V523 = "J5A(2):CASH FRM LOAN/SAV"                                           
   V524 = "J6(2):BANKCARD/STOR ACCT"                                           
   V525 = "J7(2):AMOUNT FINANCED"                                              
   V526 = "J8(2):MAKING PAYMENT NOW"                                           
   V527 = "J9(2):PAYMENT AMOUNT"                                               
   V528 = "J9(2):PAYMENT FREQUENCY"                                            
   V529 = "J10(2):CRDT INSURN W/PMT"                                           
   V530 = "J11(2):#PMTS AGREED ON"                                             
   V531 = "J12(2):PAY BNK/FC/CU/STR"                                           
   V532 = "J13(2):CREDIT THRU STORE"                                           
   V533 = "J14(2):BUY OTHR HH DRBLS"                                           
   V534 = "J2(3):ITEM CODE"                                                    
   V535 = "J3(3):WHEN BOUGHT-MONTH"                                            
   V536 = "J3(3):WHEN BOUGHT-YEAR"                                             
   V537 = "J4(3):COST OF ITEM"                                                 
   V538 = "J5(3):CASH/CRDT/CHG ACCT"                                           
   V539 = "J5A(3):CASH FRM LOAN/SAV"                                           
   V540 = "J6(3):BANKCARD/STOR ACCT"                                           
   V541 = "J7(3):AMOUNT FINANCED"                                              
   V542 = "J8(3):MAKING PAYMENT NOW"                                           
   V543 = "J9(3):PAYMENT AMOUNT"                                               
   V544 = "J9(3):PAYMENT FREQUENCY"                                            
   V545 = "J10(3):CRDT INSURN W/PMT"                                           
   V546 = "J11(3):#PMTS AGREED ON"                                             
   V547 = "J12(3):PAY BNK/FC/CU/STR"                                           
   V548 = "J13(3):CREDIT THRU STORE"                                           
   V601 = "K1:HOBBY AND REC ITEMS-76/77"                                       
   V602 = "K2(1):ITEM CODE"                                                    
   V603 = "K3(1):WHEN BOUGHT-MONTH"                                            
   V604 = "K3(1):WHEN BOUGHT-YEAR"                                             
   V605 = "K4(1):COST OF ITEM"                                                 
   V606 = "K5(1):CASH/CRDT/CHG ACCT"                                           
   V607 = "K5A(1):CASH FRM LOAN/SAV"                                           
   V608 = "K6(1):BANKCARD/STOR ACCT"                                           
   V609 = "K7(1):AMOUNT FINANCED"                                              
   V610 = "K8(1):MAKING PAYMENT NOW"                                           
   V611 = "K9(1):PAYMENT AMOUNT"                                               
   V612 = "K9(1):PAYMENT FREQUENCY"                                            
   V613 = "K10(1):CRDT INSURN W/PMT"                                           
   V614 = "K11(1):#PMTS AGREED ON"                                             
   V615 = "K12(1):PAY BNK/FC/CU/STR"                                           
   V616 = "K13(1):CREDIT THRU STORE"                                           
   V617 = "K14(1):BUY OTHER H AND R"                                           
   V618 = "K2(2):ITEM CODE"                                                    
   V619 = "K3(2):WHEN BOUGHT-MONTH"                                            
   V620 = "K3(2):WHEN BOUGHT-YEAR"                                             
   V621 = "K4(2):COST OF ITEM"                                                 
   V622 = "K5(2):CASH/CRDT/CHG ACCT"                                           
   V623 = "K5A(2):CASH FRM LOAN/SAV"                                           
   V624 = "K6(2):BANKCARD/STOR ACCT"                                           
   V625 = "K7(2):AMOUNT FINANCED"                                              
   V626 = "K8(2):MAKING PAYMENT NOW"                                           
   V627 = "K9(2):PAYMENT AMOUNT"                                               
   V628 = "K9(2):PAYMENT FREQUENCY"                                            
   V629 = "K10(2):CRDT INSURN W/PMT"                                           
   V630 = "K11(2):#PMTS AGREED ON"                                             
   V631 = "K12(2):PAY BNK/FC/CU/STR"                                           
   V632 = "K13(2):CREDIT THRU STORE"                                           
   V633 = "K14(2):BUY OTHER H AND R"                                           
   V634 = "K2(3):ITEM CODE"                                                    
   V635 = "K3(3):WHEN BOUGHT-MONTH"                                            
   V636 = "K3(3):WHEN BOUGHT-YEAR"                                             
   V637 = "K4(3):COST OF ITEM"                                                 
   V638 = "K5(3):CASH/CRDT/CHG ACCT"                                           
   V639 = "K5A(3):CASH FRM LOAN/SAV"                                           
   V640 = "K6(3):BANKCARD/STOR ACCT"                                           
   V641 = "K7(3):AMOUNT FINANCED"                                              
   V642 = "K8(3):MAKING PAYMENT NOW"                                           
   V643 = "K9(3):PAYMENT AMOUNT"                                               
   V644 = "K9(3):PAYMENT FREQUENCY"                                            
   V645 = "K10(3):CRDT INSURN W/PMT"                                           
   V646 = "K11(3):#PMTS AGREED ON"                                             
   V647 = "K12(3):PAY BNK/FC/CU/STR"                                           
   V648 = "K13(3):CREDIT THRU STORE"                                           
   V701 = "L1:R AND FU HV OTHR PMT/DBTS"                                       
   V702 = "L2(1):ITEM CODE"                                                    
   V703 = "L3(1):PAYMENT AMOUNT"                                               
   V704 = "L3(1):PAYMENT FREQUENCY"                                            
   V705 = "L4(1):ORIG AMOUNT OWED"                                             
   V706 = "L5(1):WHN TK ON DEBT-MTH"                                           
   V707 = "L5(1):WHN TK ON DEBT-YR"                                            
   V708 = "L6(1):#PMTS AGREED ON"                                              
   V709 = "L7(1):CRDT INSURN W/PMT"                                            
   V710 = "L8(1):PAY BNK/FC/CU/STR"                                            
   V711 = "L9(1):WHO ARRANGED CRDT"                                            
   V712 = "L10(1):REG PMT-TRAVL MED"                                           
   V713 = "L11(1):OTHER REG PMTS"                                              
   V714 = "L2(2):ITEM CODE"                                                    
   V715 = "L3(2):PAYMENT AMOUNT"                                               
   V716 = "L3(2):PAYMENT FREQUENCY"                                            
   V717 = "L4(2):ORIG AMOUNT OWED"                                             
   V718 = "L5(2):WHN TK ON DEBT-MTH"                                           
   V719 = "L5(2):WHN TK ON DEBT-YR"                                            
   V720 = "L6(2):#PMTS AGREED ON"                                              
   V721 = "L7(2):CRDT INSURN W/PMT"                                            
   V722 = "L8(2):PAY BNK/FC/CU/STR"                                            
   V723 = "L9(2):WHO ARRANGED CRDT"                                            
   V724 = "L10(2):REG PMT-TRAVL MED"                                           
   V725 = "L11(2):OTHER REG PMTS"                                              
   V726 = "L2(3):ITEM CODE"                                                    
   V727 = "L3(3):PAYMENT AMOUNT"                                               
   V728 = "L3(3):PAYMENT FREQUENCY"                                            
   V729 = "L4(3):ORIG AMOUNT OWED"                                             
   V730 = "L5(3):WHN TK ON DEBT-MTH"                                           
   V731 = "L5(3):WHN TK ON DEBT-YR"                                            
   V732 = "L6(3):#PMTS AGREED ON"                                              
   V733 = "L7(3):CRDT INSURN W/PMT"                                            
   V734 = "L8(3):PAY BNK/FC/CU/STR"                                            
   V735 = "L9(3):WHO ARRANGED CRDT"                                            
   V736 = "L12:DEBT-IRREGULAR PMTS"                                            
   V737 = "L13:ITEM CODE"                                                      
   V738 = "L14(1):PAY BNK/INS/INDIV"                                           
   V739 = "L15(1):ORIG AMOUNT OWED"                                            
   V740 = "L16(1):ANNUAL %RATE INTR"                                           
   V741 = "L17(1):OTHER IRREG PMTS"                                            
   V742 = "L13(2):ITEM CODE"                                                   
   V743 = "L14(2):PAY BNK/INS/INDIV"                                           
   V744 = "L15(2):ORIG AMOUNT OWED"                                            
   V745 = "L16(2):ANNUAL %RATE INTR"                                           
   V746 = "L17(2):OTHER IRREG PMTS"                                            
   V747 = "L13(3):ITEM CODE"                                                   
   V748 = "L14(3):PAY BNK/INS/INDIV"                                           
   V749 = "L15(3):ORIG AMOUNT OWED"                                            
   V750 = "L16(3):ANNUAL %RATE INTR"                                           
   V751 = "L18:CHKPT-HOME OWNERSHIP"                                           
   V752 = "L19:HOME=LOAN SECURITY?"                                            
   V753 = "L19A(1):LOAN CODE"                                                  
   V754 = "L19A(2):LOAN CODE"                                                  
   V801 = "N1:CHKPT-OUTSTNDG LOANS"                                            
   V802 = "N2:MOST RECENT LOAN CODE"                                           
   V803 = "N3:ANNUAL %RATE  INTERST"                                           
   V804 = "N3A:CHKBOX-CHECK RECORDS"                                           
   V805 = "N4:NTRST INFO ORAL/WRTN"                                            
   V806 = "N5:WHN FIND OUT NTRST RT"                                           
   V807 = "N6:TOT $ AMT INT/FIN CHG"                                           
   V808 = "N6A:CHKBOX-CHECK RECORDS"                                           
   V809 = "N6B:WHO CALC INRST RATE"                                            
   V810 = "N7:HOW SATISFIED W/LOAN"                                            
   V811 = "N7A(1):WHY SATIS W/LOAN"                                            
   V812 = "N7A(2):WHY SATIS W/LOAN"                                            
   V813 = "N8(1):WHY THIS TYPE INST"                                           
   V814 = "N8(2):WHY THIS TYPE INST"                                           
   V815 = "N8A(1):WHY PRTCLR INSTN"                                            
   V816 = "N8A(2):WHY PRTCLR INSTN"                                            
   V817 = "N9:OBTN CRDT THERE PREV"                                            
   V818 = "N9A:# TIMES OBTN CRDT"                                              
   V819 = "N10:HOW SATIS PREV CRDT"                                            
   V820 = "N11:GET INF RE:OTHR CRDT"                                           
   V821 = "N11A(1):WHAT DID YOU DO"                                            
   V822 = "N11A(2):WHAT DID YOU DO"                                            
   V823 = "N11B(1):KIND INFO WANTED"                                           
   V824 = "N11B(2):KIND INFO WANTED"                                           
   V825 = "N11C:ABLE GET INF WANTED"                                           
   V826 = "N12:GET TR-IN-LND STMT"                                             
   V827 = "N13:KEEP TR-IN-LND STMT"                                            
   V828 = "N14:LOAN HAVE CRDT INSUR"                                           
   V829 = "N15:KIND CREDIT INSURANC"                                           
   V830 = "N16:INSURN THRU CREDITOR"                                           
   V831 = "N17(1):WHY INS FROM OTHR"                                           
   V832 = "N17(2):WHY INS FROM OTHR"                                           
   V833 = "N18:SEP CHRG FOR INSURN"                                            
   V834 = "N19:CRDTR MENTION INSURN"                                           
   V835 = "N20:INS-MK DIFF GET LOAN"                                           
   V836 = "N21:INS-GET BETTER TERMS"                                           
   V837 = "N22:CRDT INSUR G/B THING"                                           
   V838 = "N23:CRDT INSUR EXPENSIVE"                                           
   V839 = "R1:CHKPT-LOAN FROM G-K"                                             
   V840 = "R2(1):ITEM FROM DEALER"                                             
   V841 = "R2(2):ITEM FROM STORE"                                              
   V842 = "R2(3):ITEM FROM CON/BLDR"                                           
   V843 = "R2(4):ITEM FROM INDIVDL"                                            
   V844 = "R2(5):ITEM FROM OTHER"                                              
   V845 = "R2A:CHKPT-# OF MENTIONS"                                            
   V846 = "R3(1):WHY THIS D/STR/CON"                                           
   V847 = "R3(1):WHY THIS D/STR/CON"                                           
   V848 = "R4:CRDT DECIDE D/STR/CON"                                           
   V849 = "R4A(1):HOW DEC D/STR/CON"                                           
   V850 = "R4A(2):HOW DEC D/STR/CON"                                           
   V851 = "R5:KNEW BEFORE BUY CRDT"                                            
   V852 = "R6:COULD HAVE PAID CASH"                                            
   V853 = "R7:CHKPT-CASH FROM G-K"                                             
   V854 = "R8:MOST RECENT CASH"                                                
   V855 = "R9(1):ITEM FROM DEALER"                                             
   V856 = "R9(2):ITEM FROM STORE"                                              
   V857 = "R9(3):ITEM FROM CON/BLDR"                                           
   V858 = "R9(4):ITEM FROM INDIVDL"                                            
   V859 = "R9(7):ITEM FROM OTHER"                                              
   V860 = "R9A:CHKPT-# OF MENTIONS"                                            
   V861 = "R10(1):WHY THIS D/ST/CON"                                           
   V862 = "R10(2):WHY THIS D/ST/CON"                                           
   V863 = "R11:WHEN DECIDE BUY CASH"                                           
   V864 = "R12:GET INF RE:OTHER DSC"                                           
   V865 = "R12A(1):SOURCE OTHR INFO"                                           
   V866 = "R12A(2):SOURCE OTHR INFO"                                           
   V867 = "R12A(3):SOURCE OTHR INFO"                                           
   V868 = "R12B(1):KIND OTHER INFO"                                            
   V869 = "R12B(2):KIND OTHER INFO"                                            
   V870 = "R12B(3):KIND OTHER INFO"                                            
   V871 = "R12C:ENOUGH INFO"                                                   
   V872 = "R13:EASY GET USEFUL INFO"                                           
   V873 = "R14:CKPT-ITM FR CNTR/D.S"                                           
   V874 = "R15:OWN BRAND BEFORE"                                               
   V875 = "R16:SATIS W/BRAND BEFORE"                                           
   V876 = "R17:#PRVS TRNS W.DL/S/CN"                                           
   V877 = "R18:SATISFIED W/PREV EXP"                                           
   V901 = "S1(1):SAV ACCT-BANK"                                                
   V902 = "S1(2):SAV ACCT-S AND LOAN"                                          
   V903 = "S1(3):SAV ACCT-MUT SV BK"                                           
   V904 = "S1(4):SAV ACCT-CREDIT UN"                                           
   V905 = "S1A:# SAVINGS ACCOUNTS"                                             
   V906 = "S1B:TOT AMT IN SAV ACCTS"                                           
   V907 = "S2:MAIN SAVINGS ACCT"                                               
   V908 = "S3(1):WHY ACCT THIS INST"                                           
   V909 = "S3(2):WHY ACCT THIS INST"                                           
   V910 = "S3A:NTRST AD INFL CHOICE"                                           
   V911 = "S4:REGULAR/CERTIFCT SVNG"                                           
   V912 = "S5:NOW ACCOUNT"                                                     
   V913 = "S6:HAVE CHECKING ACCNTS"                                            
   V914 = "S6A:# CHECKING ACCTS"                                               
   V915 = "S6B:AMT IN CHECK ACCTS"                                             
   V916 = "S6C:AFFECT BY WDRAW-WAIT"                                           
   V917 = "S7A:OWN SAVINGS BONDS"                                              
   V918 = "S8A:VALUE SAVINGS BONDS"                                            
   V919 = "S7B:OWN OTHER FED BONDS"                                            
   V920 = "S8B:VALUE OTHER FED BOND"                                           
   V921 = "S7C:OWN STATE/MUNI BONDS"                                           
   V922 = "S8C:VALUE ST/MUNI BONDS"                                            
   V923 = "S7D:OWN CORP BONDS"                                                 
   V924 = "S8D:VALUE CORP BONDS"                                               
   V925 = "S7E:OWN CORP STOCK"                                                 
   V926 = "S8E:VALUE CORP STOCK"                                               
   V927 = "S7F:OWN MUTUAL FUNDS"                                               
   V928 = "S8F:VALUE MUTUAL FUNDS"                                             
   V929 = "S7G:OWN STOCK-INVST CLUB"                                           
   V930 = "S8G:VALUE STOCK-INVST CL"                                           
   V931 = "S7H:OWN CERT OF DEPOSIT"                                            
   V932 = "S8H:VALUE CERT OF DPST"                                             
   V933 = "S9:OWN REAL ESTATE"                                                 
   V934 = "S9A(1):KIND REAL ESTATE"                                            
   V935 = "S9B(1):VALUE REAL ESTATE"                                           
   V936 = "S9C(1):AMT OWED ON PRPTY"                                           
   V937 = "S9A(2):KIND REAL ESTATE"                                            
   V938 = "S9B(2):VALUE REAL ESTATE"                                           
   V939 = "S9C(2):AMT OWED ON PRPTY"                                           
   V940 = "S9A(3):KIND REAL ESTATE"                                            
   V941 = "S9B(3):VALUE REAL ESTATE"                                           
   V942 = "S9C(3):AMT OWED ON PRPTY"                                           
   V943 = "S10:MEMBER CREDIT UNION"                                            
   V944 = "S11:NET SAVINGS LAST YR"                                            
   V945 = "S11A:LRG/SM NET SAVINGS"                                            
   V946 = "S12:WAGE INCOME OF R"                                               
   V947 = "S13-S14A:TOT FAM INCOME"                                            
   V948 = "HOUSEHOLD COMPOSITION"                                              
   V949 = "RELATION OF R TO FU HEAD"                                           
   V950 = "MARITAL STATUS OF FU HD"                                            
   V951 = "SEX OF FU HEAD"                                                     
   V952 = "AGE OF FU HEAD"                                                     
   V953 = "AGE OF WIFE OF FU HEAD"                                             
   V954 = "#PERSONS 18 OR OLDER-FU"                                            
   V955 = "ADULTS 65+ IN FU"                                                   
   V956 = "#PERSONS 17 OR YOUNGR-FU"                                           
   V957 = "AGE YOUNGEST CHILD<19-FU"                                           
   V958 = "AGE OLDEST CHILD<19-FU"                                             
   V959 = "X1:SEX OF R"                                                        
   V960 = "X2:RELTN OF R TO HH"                                                
   V961 = "X3:RACIAL/ETHNIC-5 GROUP"                                           
   V962 = "X4:R'S UNDERSTANDING"                                               
   V963 = "X5:TYPE OF STRUCTURE"                                               
   V964 = "X6:R (IN)ARTICULATE"                                                
   V965 = "X7:R SUSPICIOUS-BEFORE"                                             
   V966 = "X7A:R SUSPICIOUS-AFTER"                                             
   V967 = "X8:R'S INTEREST IN INTVW"                                           
   V968 = "X9:DID R HURRY ANSWERS"                                             
   V969 = "X10:OTHERS PRESENT INTVW"                                           
   V970 = "X11:R RFR TO LOAN PAPERS"                                           
   V1001 = "ITEM CODE"                                                         
   V1002 = "MONTH OF PURCHASE"                                                 
   V1003 = "YEAR OF PURCHASE"                                                  
   V1004 = "COST OF ITEM"                                                      
   V1005 = "CASH/CREDIT/CHARGE ACCT"                                           
   V1006 = "CASH FROM LOAN/SAVINGS"                                            
   V1007 = "BANKCARD/STORE ACCOUNT"                                            
   V1008 = "AMOUNT FINANCED"                                                   
   V1009 = "MAKING PAYMENTS NOW"                                               
   V1010 = "PAYMENT AMOUNT"                                                    
   V1011 = "PAYMENT FREQUENCY"                                                 
   V1012 = "CREDIT INSUR INCL W/PMT"                                           
   V1013 = "# PAYMENTS AGREED ON"                                              
   V1014 = "PY BNK/FNC.C/CU/DL/ST/MD"                                          
   V1015 = "CREDIT ARRANGEMENTS"                                               
   V1021 = "TRANSACTION IN N"                                                  
   V1022 = "N3-RECODED RATE"                                                   
   V1023 = "RATE AWARENESS-NCCF"                                               
   V1024 = "HAVE BANK CARD?"                                                   
   V1025 = "HAVE STORE CARD?"                                                  
   V1026 = "APR-CARDS"                                                         
   V1027 = "AGE OF R"                                                          
   V1028 = "INSTITUTION IN N"                                                  
   V1029 = "RATES-HYPOTHETICAL CASE"                                           
   V1030 = "RECODE-HYPO.RATES"                                                 
   V1031 = "APR FROM EST CHG"                                                  
   V1032 = "APR FROM 1/2 EST CHG"                                              
   V1033 = "BRACKET V1031,APR"                                                 
   V1034 = "BRACKET V1032, APR/2"                                              
   V1035 = "BRACKET V1033"                                                     
   V1036 = "BRACKET V1034";                                                    
                                                                               
* USER-DEFINED MISSING VALUE RECODE TO SAS SYSMIS;                             
                                                                               
/*                                                                             
   IF V3 GE 0000008 OR V3=0000009 THEN V3=.;                                   
   IF V4 GE 0000098 OR V4=0000000 THEN V4=.;                                   
   IF V5 GE 0000098 OR V5=0000000 THEN V5=.;                                   
   IF V6 GE 0000008 OR V6=0000009 THEN V6=.;                                   
   IF V7 GE 0000008 OR V7=0000009 THEN V7=.;                                   
   IF V8 GE 0000008 OR V8=0000009 THEN V8=.;                                   
   IF V9 GE 0000008 OR V9=0000009 THEN V9=.;                                   
   IF V10 GE 0000008 OR V10=0000009 THEN V10=.;                                
   IF V11 GE 0000008 OR V11=0000009 THEN V11=.;                                
   IF V12 GE 0000008 OR V12=0000009 THEN V12=.;                                
   IF V13 GE 0000008 OR V13=0000009 THEN V13=.;                                
   IF V14 GE 0000008 OR V14=0000009 THEN V14=.;                                
   IF V15 GE 0000008 OR V15=0000009 THEN V15=.;                                
   IF V16 GE 0000008 OR V16=0000009 THEN V16=.;                                
   IF V17 GE 0000008 OR V17=0000009 THEN V17=.;                                
   IF V18 GE 0000008 OR V18=0000009 THEN V18=.;                                
   IF V19 GE 0000008 OR V19=0000009 THEN V19=.;                                
   IF V20 GE 0000008 OR V20=0000009 THEN V20=.;                                
   IF V21 GE 0000098 OR V21=0000000 THEN V21=.;                                
   IF V22 GE 0000098 OR V22=0000000 THEN V22=.;                                
   IF V23 GE 0000098 OR V23=0000000 THEN V23=.;                                
   IF V24 GE 0000098 OR V24=0000000 THEN V24=.;                                
   IF V25 GE 0000098 OR V25=0000000 THEN V25=.;                                
   IF V26 GE 0000008 OR V26=0000009 THEN V26=.;                                
   IF V27 GE 0000008 OR V27=0000009 THEN V27=.;                                
   IF V28 GE 0000008 OR V28=0000009 THEN V28=.;                                
   IF V29 GE 0000008 OR V29=0000009 THEN V29=.;                                
   IF V30 GE 0000008 OR V30=0000009 THEN V30=.;                                
   IF V31 GE 0000008 OR V31=0000009 THEN V31=.;                                
   IF V32 GE 0000008 OR V32=0000009 THEN V32=.;                                
   IF V33 GE 0000008 OR V33=0000009 THEN V33=.;                                
   IF V34 GE 0000008 OR V34=0000009 THEN V34=.;                                
   IF V35 GE 0000008 OR V35=0000009 THEN V35=.;                                
   IF V36 GE 0000008 OR V36=0000009 THEN V36=.;                                
   IF V37 GE 0000008 OR V37=0000009 THEN V37=.;                                
   IF V38 GE 0000008 OR V38=0000009 THEN V38=.;                                
   IF V39 GE 0000008 OR V39=0000009 THEN V39=.;                                
   IF V40 GE 0000008 OR V40=0000009 THEN V40=.;                                
   IF V41 GE 0000008 OR V41=0000009 THEN V41=.;                                
   IF V42 GE 0000008 OR V42=0000009 THEN V42=.;                                
   IF V43 GE 0000008 OR V43=0000009 THEN V43=.;                                
   IF V44 GE 0000008 OR V44=0000009 THEN V44=.;                                
   IF V45 GE 0000008 OR V45=0000009 THEN V45=.;                                
   IF V46 GE 0000008 OR V46=0000009 THEN V46=.;                                
   IF V47 GE 0000008 OR V47=0000009 THEN V47=.;                                
   IF V48 GE 0000008 OR V48=0000009 THEN V48=.;                                
   IF V49 GE 0000008 OR V49=0000009 THEN V49=.;                                
   IF V50 GE 0000008 OR V50=0000009 THEN V50=.;                                
   IF V51 GE 0000098 OR V51=0000000 THEN V51=.;                                
   IF V52 GE 0000098 OR V52=0000000 THEN V52=.;                                
   IF V53 GE 0000098 OR V53=0000000 THEN V53=.;                                
   IF V54 GE 0000008 OR V54=0000009 THEN V54=.;                                
   IF V55 GE 0000008 OR V55=0000000 THEN V55=.;                                
   IF V56 GE 0000008 OR V56=0000000 THEN V56=.;                                
   IF V57 GE 0000008 OR V57=0000000 THEN V57=.;                                
   IF V58 GE 0000008 OR V58=0000000 THEN V58=.;                                
   IF V59 GE 0000098 OR V59=0000000 THEN V59=.;                                
   IF V60 GE 0000098 OR V60=0000000 THEN V60=.;                                
   IF V61 GE 0000008 OR V61=0000009 THEN V61=.;                                
   IF V62 GE 0000008 OR V62=0000009 THEN V62=.;                                
   IF V63 GE 0000008 OR V63=0000009 THEN V63=.;                                
   IF V64 GE 0000008 OR V64=0000009 THEN V64=.;                                
   IF V65 GE 0000008 OR V65=0000009 THEN V65=.;                                
   IF V66 GE 0000008 OR V66=0000000 THEN V66=.;                                
   IF V67 GE 0000008 OR V67=0000000 THEN V67=.;                                
   IF V68 GE 0000008 OR V68=0000000 THEN V68=.;                                
   IF V69 GE 0000008 OR V69=0000000 THEN V69=.;                                
   IF V70 GE 0000098 OR V70=0000000 THEN V70=.;                                
   IF V71 GE 0000098 OR V71=0000000 THEN V71=.;                                
   IF V72 GE 0009998 OR V72=0009999 THEN V72=.;                                
   IF V73 GE 0000009 OR V73=0000009 THEN V73=.;                                
   IF V74 GE 0000098 OR V74=0000099 THEN V74=.;                                
   IF V75 GE 0000009 OR V75=0000009 THEN V75=.;                                
   IF V76 GE 0000008 OR V76=0000009 THEN V76=.;                                
   IF V77 GE 0000008 OR V77=0000009 THEN V77=.;                                
   IF V78 GE 0000008 OR V78=0000009 THEN V78=.;                                
   IF V79 GE 0000008 OR V79=0000009 THEN V79=.;                                
   IF V80 GE 0000008 OR V80=0000009 THEN V80=.;                                
   IF V81 GE 0000098 OR V81=0000000 THEN V81=.;                                
   IF V82 GE 0000098 OR V82=0000000 THEN V82=.;                                
   IF V83 GE 0000008 OR V83=0000009 THEN V83=.;                                
   IF V84 GE 0000008 OR V84=0000009 THEN V84=.;                                
   IF V85 GE 0000098 OR V85=0000000 THEN V85=.;                                
   IF V86 GE 0000098 OR V86=0000000 THEN V86=.;                                
   IF V87 GE 0000098 OR V87=0000000 THEN V87=.;                                
   IF V88 GE 0000098 OR V88=0000000 THEN V88=.;                                
   IF V89 GE 0000008 OR V89=0000009 THEN V89=.;                                
   IF V90 GE 0000098 OR V90=0000000 THEN V90=.;                                
   IF V91 GE 0000098 OR V91=0000000 THEN V91=.;                                
   IF V101 GE 0000008 OR V101=0000009 THEN V101=.;                             
   IF V102 GE 0000098 OR V102=0000000 THEN V102=.;                             
   IF V103 GE 0000098 OR V103=0000000 THEN V103=.;                             
   IF V104 GE 0000008 OR V104=0000000 THEN V104=.;                             
   IF V105 GE 0000098 OR V105=0000000 THEN V105=.;                             
   IF V106 GE 0000098 OR V106=0000000 THEN V106=.;                             
   IF V107 GE 0000008 OR V107=0000000 THEN V107=.;                             
   IF V108 GE 0000008 OR V108=0000000 THEN V108=.;                             
   IF V109 GE 0000098 OR V109=0000000 THEN V109=.;                             
   IF V110 GE 0000098 OR V110=0000000 THEN V110=.;                             
   IF V111 GE 0000008 OR V111=0000000 THEN V111=.;                             
   IF V112 GE 0000098 OR V112=0000000 THEN V112=.;                             
   IF V113 GE 0000098 OR V113=0000000 THEN V113=.;                             
   IF V114 GE 0000008 OR V114=0000000 THEN V114=.;                             
   IF V115 GE 0000008 OR V115=0000009 THEN V115=.;                             
   IF V116 GE 0000098 OR V116=0000000 THEN V116=.;                             
   IF V117 GE 0000098 OR V117=0000000 THEN V117=.;                             
   IF V118 GE 0000008 OR V118=0000000 THEN V118=.;                             
   IF V119 GE 0000098 OR V119=0000000 THEN V119=.;                             
   IF V120 GE 0000098 OR V120=0000000 THEN V120=.;                             
   IF V121 GE 0000008 OR V121=0000000 THEN V121=.;                             
   IF V122 GE 0000098 OR V122=0000000 THEN V122=.;                             
   IF V123 GE 0000098 OR V123=0000000 THEN V123=.;                             
   IF V124 GE 0000008 OR V124=0000000 THEN V124=.;                             
   IF V125 GE 0000098 OR V125=0000000 THEN V125=.;                             
   IF V126 GE 0000098 OR V126=0000000 THEN V126=.;                             
   IF V127 GE 0000008 OR V127=0000009 THEN V127=.;                             
   IF V128 GE 0000098 OR V128=0000000 THEN V128=.;                             
   IF V129 GE 0000098 OR V129=0000000 THEN V129=.;                             
   IF V130 GE 0000098 OR V130=0000000 THEN V130=.;                             
   IF V131 GE 0000098 OR V131=0000000 THEN V131=.;                             
   IF V132 GE 0000008 OR V132=0000000 THEN V132=.;                             
   IF V133 GE 0000098 OR V133=0000000 THEN V133=.;                             
   IF V134 GE 0000098 OR V134=0000000 THEN V134=.;                             
   IF V135 GE 0000009 OR V135=0000009 THEN V135=.;                             
   IF V136 GE 0000008 OR V136=0000009 THEN V136=.;                             
   IF V137 GE 0000098 OR V137=0000000 THEN V137=.;                             
   IF V138 GE 0000098 OR V138=0000000 THEN V138=.;                             
   IF V139 GE 0000008 OR V139=0000009 THEN V139=.;                             
   IF V140 GE 0000008 OR V140=0000009 THEN V140=.;                             
   IF V141 GE 0000008 OR V141=0000009 THEN V141=.;                             
   IF V142 GE 0000008 OR V142=0000009 THEN V142=.;                             
   IF V143 GE 0000098 OR V143=0000000 THEN V143=.;                             
   IF V144 GE 0000098 OR V144=0000000 THEN V144=.;                             
   IF V145 GE 0000098 OR V145=0000000 THEN V145=.;                             
   IF V146 GE 0000098 OR V146=0000000 THEN V146=.;                             
   IF V147 GE 0000008 OR V147=0000009 THEN V147=.;                             
   IF V148 GE 0000008 OR V148=0000009 THEN V148=.;                             
   IF V149 GE 0000008 OR V149=0000009 THEN V149=.;                             
   IF V150 GE 0000008 OR V150=0000009 THEN V150=.;                             
   IF V151 GE 0000008 OR V151=0000009 THEN V151=.;                             
   IF V152 GE 0000008 OR V152=0000009 THEN V152=.;                             
   IF V153 GE 0000008 OR V153=0000009 THEN V153=.;                             
   IF V154 GE 0000008 OR V154=0000009 THEN V154=.;                             
   IF V155 GE 0000008 OR V155=0000009 THEN V155=.;                             
   IF V156 GE 0000008 OR V156=0000009 THEN V156=.;                             
   IF V157 GE 0000008 OR V157=0000009 THEN V157=.;                             
   IF V158 GE 0000008 OR V158=0000009 THEN V158=.;                             
   IF V159 GE 0000008 OR V159=0000009 THEN V159=.;                             
   IF V160 GE 0000008 OR V160=0000009 THEN V160=.;                             
   IF V161 GE 0000008 OR V161=0000009 THEN V161=.;                             
   IF V162 GE 0000008 OR V162=0000009 THEN V162=.;                             
   IF V163 GE 0000008 OR V163=0000009 THEN V163=.;                             
   IF V164 GE 0000008 OR V164=0000009 THEN V164=.;                             
   IF V165 GE 0000008 OR V165=0000009 THEN V165=.;                             
   IF V166 GE 0000008 OR V166=0000009 THEN V166=.;                             
   IF V167 GE 0000008 OR V167=0000009 THEN V167=.;                             
   IF V168 GE 0000008 OR V168=0000009 THEN V168=.;                             
   IF V169 GE 0000008 OR V169=0000009 THEN V169=.;                             
   IF V170 GE 0000008 OR V170=0000009 THEN V170=.;                             
   IF V171 GE 0000008 OR V171=0000009 THEN V171=.;                             
   IF V172 GE 0000008 OR V172=0000009 THEN V172=.;                             
   IF V173 GE 0000008 OR V173=0000009 THEN V173=.;                             
   IF V174 GE 0000008 OR V174=0000009 THEN V174=.;                             
   IF V175 GE 0000008 OR V175=0000009 THEN V175=.;                             
   IF V176 GE 0000008 OR V176=0000009 THEN V176=.;                             
   IF V177 GE 0000008 OR V177=0000009 THEN V177=.;                             
   IF V178 GE 0000008 OR V178=0000009 THEN V178=.;                             
   IF V179 GE 0000008 OR V179=0000009 THEN V179=.;                             
   IF V180 GE 0000008 OR V180=0000009 THEN V180=.;                             
   IF V181 GE 0000008 OR V181=0000009 THEN V181=.;                             
   IF V182 GE 0000008 OR V182=0000009 THEN V182=.;                             
   IF V183 GE 0000008 OR V183=0000009 THEN V183=.;                             
   IF V184 GE 0000098 OR V184=0000000 THEN V184=.;                             
   IF V185 GE 0000098 OR V185=0000000 THEN V185=.;                             
   IF V186 GE 0000008 OR V186=0000000 THEN V186=.;                             
   IF V187 GE 0000008 OR V187=0000009 THEN V187=.;                             
   IF V188 GE 0000098 OR V188=0000000 THEN V188=.;                             
   IF V189 GE 0000098 OR V189=0000000 THEN V189=.;                             
   IF V190 GE 0000098 OR V190=0000000 THEN V190=.;                             
   IF V191 GE 0000003 OR V191=0000009 THEN V191=.;                             
   IF V192 GE 0000008 OR V192=0000000 THEN V192=.;                             
   IF V193 GE 0000098 OR V193=0000000 THEN V193=.;                             
   IF V194 GE 0000098 OR V194=0000000 THEN V194=.;                             
   IF V195 GE 0000008 OR V195=0000000 THEN V195=.;                             
   IF V196 GE 0000098 OR V196=0000000 THEN V196=.;                             
   IF V197 GE 0000098 OR V197=0000000 THEN V197=.;                             
   IF V198 GE 0000008 OR V198=0000000 THEN V198=.;                             
   IF V201 GE 0000008 OR V201=0000009 THEN V201=.;                             
   IF V202 GE 0000008 OR V202=0000009 THEN V202=.;                             
   IF V203 GE 0000098 OR V203=0000000 THEN V203=.;                             
   IF V204 GE 0000008 OR V204=0000000 THEN V204=.;                             
   IF V205 GE 0009996 OR V205=0000000 THEN V205=.;                             
   IF V206 GE 0009996 OR V206=0000000 THEN V206=.;                             
   IF V207 GE 0000096 OR V207=0000000 THEN V207=.;                             
   IF V208 GE 0000008 OR V208=0000000 THEN V208=.;                             
   IF V209 GE 0009996 OR V209=0000000 THEN V209=.;                             
   IF V210 GE 0009996 OR V210=0000000 THEN V210=.;                             
   IF V211 GE 0000096 OR V211=0000000 THEN V211=.;                             
   IF V212 GE 0000008 OR V212=0000000 THEN V212=.;                             
   IF V213 GE 0009996 OR V213=0000000 THEN V213=.;                             
   IF V214 GE 0009996 OR V214=0000000 THEN V214=.;                             
   IF V215 GE 0000096 OR V215=0000000 THEN V215=.;                             
   IF V216 GE 0000008 OR V216=0000000 THEN V216=.;                             
   IF V217 GE 0009996 OR V217=0000000 THEN V217=.;                             
   IF V218 GE 0009996 OR V218=0000000 THEN V218=.;                             
   IF V219 GE 0000096 OR V219=0000000 THEN V219=.;                             
   IF V220 GE 0000008 OR V220=0000000 THEN V220=.;                             
   IF V221 GE 0009996 OR V221=0000000 THEN V221=.;                             
   IF V222 GE 0009996 OR V222=0000000 THEN V222=.;                             
   IF V223 GE 0000008 OR V223=0000000 THEN V223=.;                             
   IF V224 GE 0000008 OR V224=0000000 THEN V224=.;                             
   IF V225 GE 0000096 OR V225=0000000 THEN V225=.;                             
   IF V226 GE 0000098 OR V226=0000000 THEN V226=.;                             
   IF V227 GE 0000008 OR V227=0000000 THEN V227=.;                             
   IF V228 GE 0000008 OR V228=0000000 THEN V228=.;                             
   IF V229 GE 0000008 OR V229=0000000 THEN V229=.;                             
   IF V230 GE 0000098 OR V230=0000000 THEN V230=.;                             
   IF V231 GE 0000098 OR V231=0000000 THEN V231=.;                             
   IF V232 GE 0000008 OR V232=0000000 THEN V232=.;                             
   IF V233 GE 0000098 OR V233=0000000 THEN V233=.;                             
   IF V234 GE 0000098 OR V234=0000000 THEN V234=.;                             
   IF V235 GE 0000008 OR V235=0000000 THEN V235=.;                             
   IF V236 GE 0009996 OR V236=0000000 THEN V236=.;                             
   IF V237 GE 0000098 OR V237=0000000 THEN V237=.;                             
   IF V238 GE 0000098 OR V238=0000000 THEN V238=.;                             
   IF V239 GE 0000008 OR V239=0000000 THEN V239=.;                             
   IF V240 GE 0000098 OR V240=0000000 THEN V240=.;                             
   IF V241 GE 0000098 OR V241=0000000 THEN V241=.;                             
   IF V242 GE 0000008 OR V242=0000000 THEN V242=.;                             
   IF V243 GE 0000098 OR V243=0000000 THEN V243=.;                             
   IF V244 GE 0000098 OR V244=0000000 THEN V244=.;                             
   IF V245 GE 0000008 OR V245=0000000 THEN V245=.;                             
   IF V246 GE 0009996 OR V246=0000000 THEN V246=.;                             
   IF V247 GE 0000098 OR V247=0000000 THEN V247=.;                             
   IF V248 GE 0000098 OR V248=0000000 THEN V248=.;                             
   IF V249 GE 0000008 OR V249=0000000 THEN V249=.;                             
   IF V250 GE 0000098 OR V250=0000000 THEN V250=.;                             
   IF V251 GE 0000098 OR V251=0000000 THEN V251=.;                             
   IF V252 GE 0000008 OR V252=0000009 THEN V252=.;                             
   IF V253 GE 0000098 OR V253=0000000 THEN V253=.;                             
   IF V254 GE 0000098 OR V254=0000000 THEN V254=.;                             
   IF V255 GE 0000008 OR V255=0000009 THEN V255=.;                             
   IF V256 GE 0000009 OR V256=0000009 THEN V256=.;                             
   IF V257 GE 0000009 OR V257=0000000 THEN V257=.;                             
   IF V258 GE 0000003 OR V258=0000000 THEN V258=.;                             
   IF V259 GE 0000008 OR V259=0000000 THEN V259=.;                             
   IF V260 GE 0000008 OR V260=0000000 THEN V260=.;                             
   IF V261 GE 0000098 OR V261=0000000 THEN V261=.;                             
   IF V262 GE 0000098 OR V262=0000000 THEN V262=.;                             
   IF V263 GE 0000006 OR V263=0000000 THEN V263=.;                             
   IF V264 GE 0000008 OR V264=0000000 THEN V264=.;                             
   IF V265 GE 0000008 OR V265=0000000 THEN V265=.;                             
   IF V266 GE 0000098 OR V266=0000000 THEN V266=.;                             
   IF V267 GE 0000098 OR V267=0000000 THEN V267=.;                             
   IF V268 GE 0000006 OR V268=0000000 THEN V268=.;                             
   IF V275 GE 0000098 OR V275=0000099 THEN V275=.;                             
   IF V276 GE 0000098 OR V276=0000000 THEN V276=.;                             
   IF V277 GE 0000098 OR V277=0000000 THEN V277=.;                             
   IF V278 GE 0000098 OR V278=0000099 THEN V278=.;                             
   IF V279 GE 0000098 OR V279=0000099 THEN V279=.;                             
   IF V280 GE 0000098 OR V280=0000099 THEN V280=.;                             
   IF V281 GE 0009998 OR V281=0009999 THEN V281=.;                             
   IF V282 GE 0000008 OR V282=0000009 THEN V282=.;                             
   IF V283 GE 0000099 OR V283=0000000 THEN V283=.;                             
   IF V284 GE 0000098 OR V284=0000000 THEN V284=.;                             
   IF V285 GE 0000098 OR V285=0000000 THEN V285=.;                             
   IF V286 GE 0000098 OR V286=0000000 THEN V286=.;                             
   IF V287 GE 0000098 OR V287=0000000 THEN V287=.;                             
   IF V288 GE 0000098 OR V288=0000000 THEN V288=.;                             
   IF V289 GE 0009998 OR V289=0000000 THEN V289=.;                             
   IF V301 GE 0000098 OR V301=0000099 THEN V301=.;                             
   IF V302 GE 0000098 OR V302=0000099 THEN V302=.;                             
   IF V303 GE 0000008 OR V303=0000009 THEN V303=.;                             
   IF V304 GE 0000009 OR V304=0000000 THEN V304=.;                             
   IF V305 GE 0009998 OR V305=0000000 THEN V305=.;                             
   IF V306 GE 0000008 OR V306=0000000 THEN V306=.;                             
   IF V307 GE 0000008 OR V307=0000000 THEN V307=.;                             
   IF V308 GE 0999996 OR V308=0000000 THEN V308=.;                             
   IF V309 GE 0000996 OR V309=0000000 THEN V309=.;                             
   IF V310 GE 0000008 OR V310=0000000 THEN V310=.;                             
   IF V311 GE 0999996 OR V311=0000000 THEN V311=.;                             
   IF V312 GE 0000008 OR V312=0000000 THEN V312=.;                             
   IF V313 GE 0000008 OR V313=0000000 THEN V313=.;                             
   IF V314 GE 0000098 OR V314=0000000 THEN V314=.;                             
   IF V315 GE 0000008 OR V315=0000000 THEN V315=.;                             
   IF V316 GE 0000008 OR V316=0000000 THEN V316=.;                             
   IF V317 GE 0999996 OR V317=0000000 THEN V317=.;                             
   IF V318 GE 0009996 OR V318=0000000 THEN V318=.;                             
   IF V319 GE 0000008 OR V319=0000000 THEN V319=.;                             
   IF V320 GE 0000008 OR V320=0000000 THEN V320=.;                             
   IF V321 GE 0000098 OR V321=0000000 THEN V321=.;                             
   IF V322 GE 0000098 OR V322=0000000 THEN V322=.;                             
   IF V323 GE 0000098 OR V323=0000000 THEN V323=.;                             
   IF V324 GE 0009998 OR V324=0000000 THEN V324=.;                             
   IF V325 GE 0000008 OR V325=0000000 THEN V325=.;                             
   IF V326 GE 0999996 OR V326=0000000 THEN V326=.;                             
   IF V327 GE 0009996 OR V327=0000000 THEN V327=.;                             
   IF V328 GE 0000008 OR V328=0000000 THEN V328=.;                             
   IF V329 GE 0000008 OR V329=0000000 THEN V329=.;                             
   IF V330 GE 0000098 OR V330=0000000 THEN V330=.;                             
   IF V331 GE 0000098 OR V331=0000000 THEN V331=.;                             
   IF V332 GE 0000098 OR V332=0000000 THEN V332=.;                             
   IF V333 GE 0009998 OR V333=0000000 THEN V333=.;                             
   IF V334 GE 0000008 OR V334=0000000 THEN V334=.;                             
   IF V335 GE 0000008 OR V335=0000000 THEN V335=.;                             
   IF V336 GE 0000008 OR V336=0000000 THEN V336=.;                             
   IF V337 GE 0000008 OR V337=0000000 THEN V337=.;                             
   IF V338 GE 0000008 OR V338=0000000 THEN V338=.;                             
   IF V339 GE 0000098 OR V339=0000000 THEN V339=.;                             
   IF V351 GE 0000008 OR V351=0000009 THEN V351=.;                             
   IF V352 GE 0000098 OR V352=0000000 THEN V352=.;                             
   IF V353 GE 0000098 OR V353=0000000 THEN V353=.;                             
   IF V354 GE 0099998 OR V354=0000000 THEN V354=.;                             
   IF V355 GE 0000008 OR V355=0000000 THEN V355=.;                             
   IF V356 GE 0000008 OR V356=0000000 THEN V356=.;                             
   IF V357 GE 0000008 OR V357=0000000 THEN V357=.;                             
   IF V358 GE 0099998 OR V358=0000000 THEN V358=.;                             
   IF V359 GE 0000008 OR V359=0000000 THEN V359=.;                             
   IF V360 GE 0000008 OR V360=0000000 THEN V360=.;                             
   IF V361 GE 0009996 OR V361=0000000 THEN V361=.;                             
   IF V362 GE 0000008 OR V362=0000000 THEN V362=.;                             
   IF V363 GE 0000096 OR V363=0000000 THEN V363=.;                             
   IF V364 GE 0000008 OR V364=0000000 THEN V364=.;                             
   IF V365 GE 0000008 OR V365=0000000 THEN V365=.;                             
   IF V366 GE 0000008 OR V366=0000000 THEN V366=.;                             
   IF V367 GE 0000098 OR V367=0000000 THEN V367=.;                             
   IF V368 GE 0000098 OR V368=0000000 THEN V368=.;                             
   IF V369 GE 0099998 OR V369=0000000 THEN V369=.;                             
   IF V370 GE 0000008 OR V370=0000000 THEN V370=.;                             
   IF V371 GE 0000008 OR V371=0000000 THEN V371=.;                             
   IF V372 GE 0000008 OR V372=0000000 THEN V372=.;                             
   IF V373 GE 0099998 OR V373=0000000 THEN V373=.;                             
   IF V374 GE 0000008 OR V374=0000000 THEN V374=.;                             
   IF V375 GE 0000008 OR V375=0000000 THEN V375=.;                             
   IF V376 GE 0009996 OR V376=0000000 THEN V376=.;                             
   IF V377 GE 0000008 OR V377=0000000 THEN V377=.;                             
   IF V378 GE 0000096 OR V378=0000000 THEN V378=.;                             
   IF V379 GE 0000008 OR V379=0000000 THEN V379=.;                             
   IF V380 GE 0000008 OR V380=0000000 THEN V380=.;                             
   IF V381 GE 0000008 OR V381=0000000 THEN V381=.;                             
   IF V382 GE 0000098 OR V382=0000000 THEN V382=.;                             
   IF V383 GE 0000098 OR V383=0000000 THEN V383=.;                             
   IF V384 GE 0099998 OR V384=0000000 THEN V384=.;                             
   IF V385 GE 0000008 OR V385=0000000 THEN V385=.;                             
   IF V386 GE 0000008 OR V386=0000000 THEN V386=.;                             
   IF V387 GE 0000008 OR V387=0000000 THEN V387=.;                             
   IF V388 GE 0099998 OR V388=0000000 THEN V388=.;                             
   IF V389 GE 0000008 OR V389=0000000 THEN V389=.;                             
   IF V390 GE 0000008 OR V390=0000000 THEN V390=.;                             
   IF V391 GE 0009996 OR V391=0000000 THEN V391=.;                             
   IF V392 GE 0000008 OR V392=0000000 THEN V392=.;                             
   IF V393 GE 0000096 OR V393=0000000 THEN V393=.;                             
   IF V394 GE 0000008 OR V394=0000000 THEN V394=.;                             
   IF V395 GE 0000008 OR V395=0000000 THEN V395=.;                             
   IF V425 GE 0000008 OR V425=0000009 THEN V425=.;                             
   IF V426 GE 0000008 OR V426=0000000 THEN V426=.;                             
   IF V427 GE 0000098 OR V427=0000000 THEN V427=.;                             
   IF V428 GE 0000098 OR V428=0000000 THEN V428=.;                             
   IF V429 GE 0000098 OR V429=0000000 THEN V429=.;                             
   IF V430 GE 0000008 OR V430=0000000 THEN V430=.;                             
   IF V431 GE 0000998 OR V431=0000000 THEN V431=.;                             
   IF V432 GE 0099996 OR V432=0000000 THEN V432=.;                             
   IF V433 GE 0000008 OR V433=0000000 THEN V433=.;                             
   IF V434 GE 0000008 OR V434=0000000 THEN V434=.;                             
   IF V435 GE 0000008 OR V435=0000000 THEN V435=.;                             
   IF V436 GE 0000008 OR V436=0000000 THEN V436=.;                             
   IF V437 GE 0000008 OR V437=0000000 THEN V437=.;                             
   IF V438 GE 0099998 OR V438=0000000 THEN V438=.;                             
   IF V439 GE 0009996 OR V439=0000000 THEN V439=.;                             
   IF V440 GE 0000008 OR V440=0000000 THEN V440=.;                             
   IF V441 GE 0000008 OR V441=0000000 THEN V441=.;                             
   IF V442 GE 0000096 OR V442=0000000 THEN V442=.;                             
   IF V443 GE 0000008 OR V443=0000000 THEN V443=.;                             
   IF V444 GE 0009996 OR V444=0000000 THEN V444=.;                             
   IF V445 GE 0000008 OR V445=0000000 THEN V445=.;                             
   IF V446 GE 0000008 OR V446=0000000 THEN V446=.;                             
   IF V447 GE 0000008 OR V447=0000000 THEN V447=.;                             
   IF V448 GE 0000098 OR V448=0000000 THEN V448=.;                             
   IF V449 GE 0000098 OR V449=0000000 THEN V449=.;                             
   IF V450 GE 0000098 OR V450=0000000 THEN V450=.;                             
   IF V451 GE 0000008 OR V451=0000000 THEN V451=.;                             
   IF V452 GE 0000998 OR V452=0000000 THEN V452=.;                             
   IF V453 GE 0099996 OR V453=0000000 THEN V453=.;                             
   IF V454 GE 0000008 OR V454=0000000 THEN V454=.;                             
   IF V455 GE 0000008 OR V455=0000000 THEN V455=.;                             
   IF V456 GE 0000008 OR V456=0000000 THEN V456=.;                             
   IF V457 GE 0000008 OR V457=0000000 THEN V457=.;                             
   IF V458 GE 0000008 OR V458=0000000 THEN V458=.;                             
   IF V459 GE 0099998 OR V459=0000000 THEN V459=.;                             
   IF V460 GE 0009996 OR V460=0000000 THEN V460=.;                             
   IF V461 GE 0000008 OR V461=0000000 THEN V461=.;                             
   IF V462 GE 0000008 OR V462=0000000 THEN V462=.;                             
   IF V463 GE 0000096 OR V463=0000000 THEN V463=.;                             
   IF V464 GE 0000008 OR V464=0000000 THEN V464=.;                             
   IF V465 GE 0009996 OR V465=0000000 THEN V465=.;                             
   IF V466 GE 0000008 OR V466=0000000 THEN V466=.;                             
   IF V467 GE 0000008 OR V467=0000000 THEN V467=.;                             
   IF V468 GE 0000008 OR V468=0000000 THEN V468=.;                             
   IF V469 GE 0000098 OR V469=0000000 THEN V469=.;                             
   IF V470 GE 0000098 OR V470=0000000 THEN V470=.;                             
   IF V471 GE 0000098 OR V471=0000000 THEN V471=.;                             
   IF V472 GE 0000008 OR V472=0000000 THEN V472=.;                             
   IF V473 GE 0000998 OR V473=0000000 THEN V473=.;                             
   IF V474 GE 0099996 OR V474=0000000 THEN V474=.;                             
   IF V475 GE 0000008 OR V475=0000000 THEN V475=.;                             
   IF V476 GE 0000008 OR V476=0000000 THEN V476=.;                             
   IF V477 GE 0000008 OR V477=0000000 THEN V477=.;                             
   IF V478 GE 0000008 OR V478=0000000 THEN V478=.;                             
   IF V479 GE 0000008 OR V479=0000000 THEN V479=.;                             
   IF V480 GE 0099998 OR V480=0000000 THEN V480=.;                             
   IF V481 GE 0009996 OR V481=0000000 THEN V481=.;                             
   IF V482 GE 0000008 OR V482=0000000 THEN V482=.;                             
   IF V483 GE 0000008 OR V483=0000000 THEN V483=.;                             
   IF V484 GE 0000096 OR V484=0000000 THEN V484=.;                             
   IF V485 GE 0000008 OR V485=0000000 THEN V485=.;                             
   IF V486 GE 0009996 OR V486=0000000 THEN V486=.;                             
   IF V487 GE 0000008 OR V487=0000000 THEN V487=.;                             
   IF V488 GE 0000008 OR V488=0000000 THEN V488=.;                             
   IF V501 GE 0000008 OR V501=0000009 THEN V501=.;                             
   IF V502 GE 0000098 OR V502=0000000 THEN V502=.;                             
   IF V503 GE 0000098 OR V503=0000000 THEN V503=.;                             
   IF V504 GE 0000098 OR V504=0000000 THEN V504=.;                             
   IF V505 GE 0099996 OR V505=0000000 THEN V505=.;                             
   IF V506 GE 0000008 OR V506=0000000 THEN V506=.;                             
   IF V507 GE 0000008 OR V507=0000000 THEN V507=.;                             
   IF V508 GE 0000008 OR V508=0000000 THEN V508=.;                             
   IF V509 GE 0099996 OR V509=0000000 THEN V509=.;                             
   IF V510 GE 0000008 OR V510=0000000 THEN V510=.;                             
   IF V511 GE 0009996 OR V511=0000000 THEN V511=.;                             
   IF V512 GE 0000008 OR V512=0000000 THEN V512=.;                             
   IF V513 GE 0000008 OR V513=0000000 THEN V513=.;                             
   IF V514 GE 0000096 OR V514=0000000 THEN V514=.;                             
   IF V515 GE 0000008 OR V515=0000000 THEN V515=.;                             
   IF V516 GE 0000008 OR V516=0000000 THEN V516=.;                             
   IF V517 GE 0000006 OR V517=0000000 THEN V517=.;                             
   IF V518 GE 0000098 OR V518=0000000 THEN V518=.;                             
   IF V519 GE 0000098 OR V519=0000000 THEN V519=.;                             
   IF V520 GE 0000098 OR V520=0000000 THEN V520=.;                             
   IF V521 GE 0099996 OR V521=0000000 THEN V521=.;                             
   IF V522 GE 0000008 OR V522=0000000 THEN V522=.;                             
   IF V523 GE 0000008 OR V523=0000000 THEN V523=.;                             
   IF V524 GE 0000008 OR V524=0000000 THEN V524=.;                             
   IF V525 GE 0099996 OR V525=0000000 THEN V525=.;                             
   IF V526 GE 0000008 OR V526=0000000 THEN V526=.;                             
   IF V527 GE 0009996 OR V527=0000000 THEN V527=.;                             
   IF V528 GE 0000008 OR V528=0000000 THEN V528=.;                             
   IF V529 GE 0000008 OR V529=0000000 THEN V529=.;                             
   IF V530 GE 0000096 OR V530=0000000 THEN V530=.;                             
   IF V531 GE 0000008 OR V531=0000000 THEN V531=.;                             
   IF V532 GE 0000008 OR V532=0000000 THEN V532=.;                             
   IF V533 GE 0000006 OR V533=0000000 THEN V533=.;                             
   IF V534 GE 0000098 OR V534=0000000 THEN V534=.;                             
   IF V535 GE 0000098 OR V535=0000000 THEN V535=.;                             
   IF V536 GE 0000098 OR V536=0000000 THEN V536=.;                             
   IF V537 GE 0099996 OR V537=0000000 THEN V537=.;                             
   IF V538 GE 0000008 OR V538=0000000 THEN V538=.;                             
   IF V539 GE 0000008 OR V539=0000000 THEN V539=.;                             
   IF V540 GE 0000008 OR V540=0000000 THEN V540=.;                             
   IF V541 GE 0099996 OR V541=0000000 THEN V541=.;                             
   IF V542 GE 0000008 OR V542=0000000 THEN V542=.;                             
   IF V543 GE 0009996 OR V543=0000000 THEN V543=.;                             
   IF V544 GE 0000008 OR V544=0000000 THEN V544=.;                             
   IF V545 GE 0000008 OR V545=0000000 THEN V545=.;                             
   IF V546 GE 0000096 OR V546=0000000 THEN V546=.;                             
   IF V547 GE 0000008 OR V547=0000000 THEN V547=.;                             
   IF V548 GE 0000008 OR V548=0000000 THEN V548=.;                             
   IF V601 GE 0000006 OR V601=0000009 THEN V601=.;                             
   IF V602 GE 0000098 OR V602=0000000 THEN V602=.;                             
   IF V603 GE 0000098 OR V603=0000000 THEN V603=.;                             
   IF V604 GE 0000098 OR V604=0000000 THEN V604=.;                             
   IF V605 GE 0099996 OR V605=0000000 THEN V605=.;                             
   IF V606 GE 0000008 OR V606=0000000 THEN V606=.;                             
   IF V607 GE 0000008 OR V607=0000000 THEN V607=.;                             
   IF V608 GE 0000008 OR V608=0000000 THEN V608=.;                             
   IF V609 GE 0099996 OR V609=0000000 THEN V609=.;                             
   IF V610 GE 0000008 OR V610=0000000 THEN V610=.;                             
   IF V611 GE 0009996 OR V611=0000000 THEN V611=.;                             
   IF V612 GE 0000008 OR V612=0000000 THEN V612=.;                             
   IF V613 GE 0000008 OR V613=0000000 THEN V613=.;                             
   IF V614 GE 0000096 OR V614=0000000 THEN V614=.;                             
   IF V615 GE 0000008 OR V615=0000000 THEN V615=.;                             
   IF V616 GE 0000008 OR V616=0000000 THEN V616=.;                             
   IF V617 GE 0000006 OR V617=0000000 THEN V617=.;                             
   IF V618 GE 0000098 OR V618=0000000 THEN V618=.;                             
   IF V619 GE 0000098 OR V619=0000000 THEN V619=.;                             
   IF V620 GE 0000098 OR V620=0000000 THEN V620=.;                             
   IF V621 GE 0099996 OR V621=0000000 THEN V621=.;                             
   IF V622 GE 0000008 OR V622=0000000 THEN V622=.;                             
   IF V623 GE 0000008 OR V623=0000000 THEN V623=.;                             
   IF V624 GE 0000008 OR V624=0000000 THEN V624=.;                             
   IF V625 GE 0099996 OR V625=0000000 THEN V625=.;                             
   IF V626 GE 0000008 OR V626=0000000 THEN V626=.;                             
   IF V627 GE 0009996 OR V627=0000000 THEN V627=.;                             
   IF V628 GE 0000008 OR V628=0000000 THEN V628=.;                             
   IF V629 GE 0000008 OR V629=0000000 THEN V629=.;                             
   IF V630 GE 0000096 OR V630=0000000 THEN V630=.;                             
   IF V631 GE 0000008 OR V631=0000000 THEN V631=.;                             
   IF V632 GE 0000008 OR V632=0000000 THEN V632=.;                             
   IF V633 GE 0000006 OR V633=0000000 THEN V633=.;                             
   IF V634 GE 0000098 OR V634=0000000 THEN V634=.;                             
   IF V635 GE 0000098 OR V635=0000000 THEN V635=.;                             
   IF V636 GE 0000098 OR V636=0000000 THEN V636=.;                             
   IF V637 GE 0099996 OR V637=0000000 THEN V637=.;                             
   IF V638 GE 0000008 OR V638=0000000 THEN V638=.;                             
   IF V639 GE 0000008 OR V639=0000000 THEN V639=.;                             
   IF V640 GE 0000008 OR V640=0000000 THEN V640=.;                             
   IF V641 GE 0099996 OR V641=0000000 THEN V641=.;                             
   IF V642 GE 0000008 OR V642=0000000 THEN V642=.;                             
   IF V643 GE 0009996 OR V643=0000000 THEN V643=.;                             
   IF V644 GE 0000008 OR V644=0000000 THEN V644=.;                             
   IF V645 GE 0000008 OR V645=0000000 THEN V645=.;                             
   IF V646 GE 0000096 OR V646=0000000 THEN V646=.;                             
   IF V647 GE 0000008 OR V647=0000000 THEN V647=.;                             
   IF V648 GE 0000008 OR V648=0000000 THEN V648=.;                             
   IF V701 GE 0000006 OR V701=0000009 THEN V701=.;                             
   IF V702 GE 0000098 OR V702=0000000 THEN V702=.;                             
   IF V703 GE 0009996 OR V703=0000000 THEN V703=.;                             
   IF V704 GE 0000008 OR V704=0000000 THEN V704=.;                             
   IF V705 GE 0099996 OR V705=0000000 THEN V705=.;                             
   IF V706 GE 0000098 OR V706=0000000 THEN V706=.;                             
   IF V707 GE 0000098 OR V707=0000000 THEN V707=.;                             
   IF V708 GE 0000096 OR V708=0000000 THEN V708=.;                             
   IF V709 GE 0000008 OR V709=0000000 THEN V709=.;                             
   IF V710 GE 0000008 OR V710=0000000 THEN V710=.;                             
   IF V711 GE 0000008 OR V711=0000000 THEN V711=.;                             
   IF V712 GE 0000006 OR V712=0000009 THEN V712=.;                             
   IF V713 GE 0000008 OR V713=0000000 THEN V713=.;                             
   IF V714 GE 0000098 OR V714=0000000 THEN V714=.;                             
   IF V715 GE 0009996 OR V715=0000000 THEN V715=.;                             
   IF V716 GE 0000008 OR V716=0000000 THEN V716=.;                             
   IF V717 GE 0099996 OR V717=0000000 THEN V717=.;                             
   IF V718 GE 0000098 OR V718=0000000 THEN V718=.;                             
   IF V719 GE 0000098 OR V719=0000000 THEN V719=.;                             
   IF V720 GE 0000096 OR V720=0000000 THEN V720=.;                             
   IF V721 GE 0000008 OR V721=0000000 THEN V721=.;                             
   IF V722 GE 0000008 OR V722=0000000 THEN V722=.;                             
   IF V723 GE 0000008 OR V723=0000000 THEN V723=.;                             
   IF V724 GE 0000008 OR V724=0000000 THEN V724=.;                             
   IF V725 GE 0000008 OR V725=0000000 THEN V725=.;                             
   IF V726 GE 0000098 OR V726=0000000 THEN V726=.;                             
   IF V727 GE 0009996 OR V727=0000000 THEN V727=.;                             
   IF V728 GE 0000008 OR V728=0000000 THEN V728=.;                             
   IF V729 GE 0099996 OR V729=0000000 THEN V729=.;                             
   IF V730 GE 0000098 OR V730=0000000 THEN V730=.;                             
   IF V731 GE 0000098 OR V731=0000000 THEN V731=.;                             
   IF V732 GE 0000096 OR V732=0000000 THEN V732=.;                             
   IF V733 GE 0000008 OR V733=0000000 THEN V733=.;                             
   IF V734 GE 0000008 OR V734=0000000 THEN V734=.;                             
   IF V735 GE 0000008 OR V735=0000000 THEN V735=.;                             
   IF V736 GE 0000006 OR V736=0000009 THEN V736=.;                             
   IF V737 GE 0000098 OR V737=0000000 THEN V737=.;                             
   IF V738 GE 0000008 OR V738=0000000 THEN V738=.;                             
   IF V739 GE 0099996 OR V739=0000000 THEN V739=.;                             
   IF V740 GE 0000096 OR V740=0000000 THEN V740=.;                             
   IF V741 GE 0000008 OR V741=0000000 THEN V741=.;                             
   IF V742 GE 0000098 OR V742=0000000 THEN V742=.;                             
   IF V743 GE 0000008 OR V743=0000000 THEN V743=.;                             
   IF V744 GE 0099996 OR V744=0000000 THEN V744=.;                             
   IF V745 GE 0000096 OR V745=0000000 THEN V745=.;                             
   IF V746 GE 0000008 OR V746=0000000 THEN V746=.;                             
   IF V747 GE 0000098 OR V747=0000000 THEN V747=.;                             
   IF V748 GE 0000008 OR V748=0000000 THEN V748=.;                             
   IF V749 GE 0099996 OR V749=0000000 THEN V749=.;                             
   IF V750 GE 0000096 OR V750=0000000 THEN V750=.;                             
   IF V751=0000009 THEN V751=.;                                                
   IF V752 GE 0000008 OR V752=0000000 THEN V752=.;                             
   IF V753 GE 0000098 OR V753=0000000 THEN V753=.;                             
   IF V754 GE 0000098 OR V754=0000000 THEN V754=.;                             
   IF V801=0000009 THEN V801=.;                                                
   IF V802 GE 0000098 OR V802=0000000 THEN V802=.;                             
   IF V803 GE 0000096 OR V803=0000000 THEN V803=.;                             
   IF V804 GE 0000008 OR V804=0000000 THEN V804=.;                             
   IF V805 GE 0000008 OR V805=0000000 THEN V805=.;                             
   IF V806 GE 0000008 OR V806=0000000 THEN V806=.;                             
   IF V807 GE 0009996 OR V807=0000000 THEN V807=.;                             
   IF V808 GE 0000009 OR V808=0000000 THEN V808=.;                             
   IF V809 GE 0000008 OR V809=0000000 THEN V809=.;                             
   IF V810 GE 0000008 OR V810=0000000 THEN V810=.;                             
   IF V811 GE 0000098 OR V811=0000000 THEN V811=.;                             
   IF V812 GE 0000098 OR V812=0000000 THEN V812=.;                             
   IF V813 GE 0000098 OR V813=0000000 THEN V813=.;                             
   IF V814 GE 0000098 OR V814=0000000 THEN V814=.;                             
   IF V815 GE 0000098 OR V815=0000000 THEN V815=.;                             
   IF V816 GE 0000098 OR V816=0000000 THEN V816=.;                             
   IF V817 GE 0000008 OR V817=0000000 THEN V817=.;                             
   IF V818 GE 0000098 OR V818=0000000 THEN V818=.;                             
   IF V819 GE 0000008 OR V819=0000000 THEN V819=.;                             
   IF V820 GE 0000008 OR V820=0000000 THEN V820=.;                             
   IF V821 GE 0000098 OR V821=0000000 THEN V821=.;                             
   IF V822 GE 0000098 OR V822=0000000 THEN V822=.;                             
   IF V823 GE 0000098 OR V823=0000000 THEN V823=.;                             
   IF V824 GE 0000098 OR V824=0000000 THEN V824=.;                             
   IF V825 GE 0000008 OR V825=0000000 THEN V825=.;                             
   IF V826 GE 0000008 OR V826=0000000 THEN V826=.;                             
   IF V827 GE 0000008 OR V827=0000000 THEN V827=.;                             
   IF V828 GE 0000008 OR V828=0000000 THEN V828=.;                             
   IF V829 GE 0000008 OR V829=0000000 THEN V829=.;                             
   IF V830 GE 0000008 OR V830=0000000 THEN V830=.;                             
   IF V831 GE 0000098 OR V831=0000000 THEN V831=.;                             
   IF V832 GE 0000098 OR V832=0000000 THEN V832=.;                             
   IF V833 GE 0000008 OR V833=0000000 THEN V833=.;                             
   IF V834 GE 0000008 OR V834=0000000 THEN V834=.;                             
   IF V835 GE 0000008 OR V835=0000000 THEN V835=.;                             
   IF V836 GE 0000008 OR V836=0000000 THEN V836=.;                             
   IF V837 GE 0000008 OR V837=0000000 THEN V837=.;                             
   IF V838 GE 0000008 OR V838=0000000 THEN V838=.;                             
   IF V839 GE 0000009 OR V839=0000000 THEN V839=.;                             
   IF V840 GE 0000008 OR V840=0000000 THEN V840=.;                             
   IF V841 GE 0000008 OR V841=0000000 THEN V841=.;                             
   IF V842 GE 0000008 OR V842=0000000 THEN V842=.;                             
   IF V843 GE 0000008 OR V843=0000000 THEN V843=.;                             
   IF V844 GE 0000008 OR V844=0000000 THEN V844=.;                             
   IF V845 GE 0000008 OR V845=0000000 THEN V845=.;                             
   IF V846 GE 0000098 OR V846=0000000 THEN V846=.;                             
   IF V847 GE 0000098 OR V847=0000000 THEN V847=.;                             
   IF V848 GE 0000008 OR V848=0000000 THEN V848=.;                             
   IF V849 GE 0000098 OR V849=0000000 THEN V849=.;                             
   IF V850 GE 0000098 OR V850=0000000 THEN V850=.;                             
   IF V851=0000000 THEN V851=.;                                                
   IF V852 GE 0000008 OR V852=0000000 THEN V852=.;                             
   IF V853=0000000 THEN V853=.;                                                
   IF V854 GE 0000098 OR V854=0000000 THEN V854=.;                             
   IF V855 GE 0000008 OR V855=0000000 THEN V855=.;                             
   IF V856 GE 0000008 OR V856=0000000 THEN V856=.;                             
   IF V857 GE 0000008 OR V857=0000000 THEN V857=.;                             
   IF V858 GE 0000008 OR V858=0000000 THEN V858=.;                             
   IF V859 GE 0000008 OR V859=0000000 THEN V859=.;                             
   IF V860 GE 0000008 OR V860=0000000 THEN V860=.;                             
   IF V861 GE 0000098 OR V861=0000000 THEN V861=.;                             
   IF V862 GE 0000098 OR V862=0000000 THEN V862=.;                             
   IF V863 GE 0000008 OR V863=0000000 THEN V863=.;                             
   IF V864 GE 0000008 OR V864=0000000 THEN V864=.;                             
   IF V865 GE 0000098 OR V865=0000000 THEN V865=.;                             
   IF V866 GE 0000098 OR V866=0000000 THEN V866=.;                             
   IF V867 GE 0000098 OR V867=0000000 THEN V867=.;                             
   IF V868 GE 0000098 OR V868=0000000 THEN V868=.;                             
   IF V869 GE 0000098 OR V869=0000000 THEN V869=.;                             
   IF V870 GE 0000098 OR V870=0000000 THEN V870=.;                             
   IF V871 GE 0000008 OR V871=0000000 THEN V871=.;                             
   IF V872 GE 0000008 OR V872=0000000 THEN V872=.;                             
   IF V873=0000000 THEN V873=.;                                                
   IF V874 GE 0000008 OR V874=0000000 THEN V874=.;                             
   IF V875 GE 0000008 OR V875=0000000 THEN V875=.;                             
   IF V876 GE 0000096 OR V876=0000000 THEN V876=.;                             
   IF V877 GE 0000008 OR V877=0000000 THEN V877=.;                             
   IF V901 GE 0000008 OR V901=0000000 THEN V901=.;                             
   IF V902 GE 0000008 OR V902=0000000 THEN V902=.;                             
   IF V903 GE 0000008 OR V903=0000000 THEN V903=.;                             
   IF V904 GE 0000008 OR V904=0000000 THEN V904=.;                             
   IF V905 GE 0000098 OR V905=0000000 THEN V905=.;                             
   IF V906 GE 0000098 OR V906=0000000 THEN V906=.;                             
   IF V907 GE 0000008 OR V907=0000000 THEN V907=.;                             
   IF V908 GE 0000098 OR V908=0000000 THEN V908=.;                             
   IF V909 GE 0000098 OR V909=0000000 THEN V909=.;                             
   IF V910 GE 0000008 OR V910=0000000 THEN V910=.;                             
   IF V911 GE 0000008 OR V911=0000000 THEN V911=.;                             
   IF V912 GE 0000008 OR V912=0000009 THEN V912=.;                             
   IF V913 GE 0000008 OR V913=0000009 THEN V913=.;                             
   IF V914 GE 0000098 OR V914=0000000 THEN V914=.;                             
   IF V915 GE 0000098 OR V915=0000000 THEN V915=.;                             
   IF V916 GE 0000008 OR V916=0000000 THEN V916=.;                             
   IF V917 GE 0000008 OR V917=0000009 THEN V917=.;                             
   IF V918 GE 0000008 OR V918=0000000 THEN V918=.;                             
   IF V919 GE 0000008 OR V919=0000009 THEN V919=.;                             
   IF V920 GE 0000098 OR V920=0000000 THEN V920=.;                             
   IF V921 GE 0000008 OR V921=0000009 THEN V921=.;                             
   IF V922 GE 0000098 OR V922=0000000 THEN V922=.;                             
   IF V923 GE 0000008 OR V923=0000009 THEN V923=.;                             
   IF V924 GE 0000098 OR V924=0000000 THEN V924=.;                             
   IF V925 GE 0000008 OR V925=0000009 THEN V925=.;                             
   IF V926 GE 0000098 OR V926=0000000 THEN V926=.;                             
   IF V927 GE 0000008 OR V927=0000009 THEN V927=.;                             
   IF V928 GE 0000098 OR V928=0000000 THEN V928=.;                             
   IF V929 GE 0000008 OR V929=0000009 THEN V929=.;                             
   IF V930 GE 0000098 OR V930=0000000 THEN V930=.;                             
   IF V931 GE 0000008 OR V931=0000009 THEN V931=.;                             
   IF V932 GE 0000098 OR V932=0000000 THEN V932=.;                             
   IF V933 GE 0000008 OR V933=0000009 THEN V933=.;                             
   IF V934 GE 0000098 OR V934=0000000 THEN V934=.;                             
   IF V935 GE 0099996 OR V935=0000000 THEN V935=.;                             
   IF V936 GE 0099996 OR V936=0000000 THEN V936=.;                             
   IF V937 GE 0000098 OR V937=0000000 THEN V937=.;                             
   IF V938 GE 0099996 OR V938=0000000 THEN V938=.;                             
   IF V939 GE 0099996 OR V939=0000000 THEN V939=.;                             
   IF V940 GE 0000098 OR V940=0000000 THEN V940=.;                             
   IF V941 GE 0099996 OR V941=0000000 THEN V941=.;                             
   IF V942 GE 0099996 OR V942=0000000 THEN V942=.;                             
   IF V943 GE 0000008 OR V943=0000009 THEN V943=.;                             
   IF V944 GE 0000008 OR V944=0000009 THEN V944=.;                             
   IF V945 GE 0000008 OR V945=0000000 THEN V945=.;                             
   IF V946 GE 0000098 OR V946=0000099 THEN V946=.;                             
   IF V947 GE 0000098 OR V947=0000099 THEN V947=.;                             
   IF V948 GE 0000098 OR V948=0000099 THEN V948=.;                             
   IF V949=0000009 THEN V949=.;        IF V950=0000009 THEN V950=.;            
   IF V951=0000009 THEN V951=.;                                                
   IF V952 GE 0000098 OR V952=0000099 THEN V952=.;                             
   IF V953 GE 0000099 OR V953=0000000 THEN V953=.;                             
   IF V954=0000099 THEN V954=.;        IF V955=0000099 THEN V955=.;            
   IF V956=0000099 THEN V956=.;                                                
   IF V957 GE 0000099 OR V957=0000000 THEN V957=.;                             
   IF V958 GE 0000099 OR V958=0000000 THEN V958=.;                             
   IF V959=0000009 THEN V959=.;        IF V960=0000009 THEN V960=.;            
   IF V961=0000009 THEN V961=.;        IF V962=0000009 THEN V962=.;            
   IF V963 GE 0000098 OR V963=0000099 THEN V963=.;                             
   IF V964=0000009 THEN V964=.;        IF V965=0000009 THEN V965=.;            
   IF V966=0000009 THEN V966=.;        IF V967=0000009 THEN V967=.;            
   IF V968=0000009 THEN V968=.;        IF V969=0000099 THEN V969=.;            
   IF V970=0000009 THEN V970=.;                                                
   IF V1001 GE 0000098 OR V1001=0000000 THEN V1001=.;                          
   IF V1002 GE 0000008 OR V1002=0000000 THEN V1002=.;                          
   IF V1003 GE 0000008 OR V1003=0000000 THEN V1003=.;                          
   IF V1004 GE 0099996 OR V1004=0000000 THEN V1004=.;                          
   IF V1005 GE 0000008 OR V1005=0000000 THEN V1005=.;                          
   IF V1006 GE 0000008 OR V1006=0000000 THEN V1006=.;                          
   IF V1007 GE 0000008 OR V1007=0000000 THEN V1007=.;                          
   IF V1008 GE 0099996 OR V1008=0000000 THEN V1008=.;                          
   IF V1009 GE 0000008 OR V1009=0000000 THEN V1009=.;                          
   IF V1010 GE 0009996 OR V1010=0000000 THEN V1010=.;                          
   IF V1011 GE 0000008 OR V1011=0000000 THEN V1011=.;                          
   IF V1012 GE 0000008 OR V1012=0000000 THEN V1012=.;                          
   IF V1013 GE 0000096 OR V1013=0000000 THEN V1013=.;                          
   IF V1014 GE 0000008 OR V1014=0000000 THEN V1014=.;                          
   IF V1015 GE 0000008 OR V1015=0000000 THEN V1015=.;                          
   IF V1021=0000000 THEN V1021=.;                                              
   IF V1022 GE 0000008 OR V1022=0000000 THEN V1022=.;                          
   IF V1023 GE 0000008 OR V1023=0000000 THEN V1023=.;                          
   IF V1024=0000009 THEN V1024=.;      IF V1025=0000009 THEN V1025=.;          
   IF V1026 GE 0000008 OR V1026=0000000 THEN V1026=.;                          
   IF V1027=0000099 THEN V1027=.;                                              
   IF V1028 GE 0000009 OR V1028=0000000 THEN V1028=.;                          
   IF V1029 GE 0000098 OR V1029=0000099 THEN V1029=.;                          
   IF V1030 GE 0000008 OR V1030=0000009 THEN V1030=.;                          
   IF V1031=09999.9 THEN V1031=.;      IF V1032=09999.9 THEN V1032=.;          
   IF V1033 GE 0000098 OR V1033=0000099 THEN V1033=.;                          
   IF V1034 GE 0000099 OR V1034=0000098 THEN V1034=.;                          
   IF V1035 GE 0000008 OR V1035=0000009 THEN V1035=.;                          
   IF V1036 GE 0000008 OR V1036=0000009 THEN V1036=.;                          
                                                                               
*/                                                                             
