	  �.  K   k820309    �
          11.1        vfV                                                                                                           
       tpcore_fvdas_mod.F90 TPCORE_FVDAS_MOD              INIT_TPCORE EXIT_TPCORE TPCORE_FVDAS FZPPM #         @                                                 
   #INIT_TPCORE%REPEAT    #INIT_TPCORE%COS    #INIT_TPCORE%SIN    #INIT_TPCORE%DBLE    #IM    #JM    #KM    #JFIRST 	   #JLAST 
   #NG    #MG    #DT    #AE    #CLAT                                      @                                 REPEAT               @                                 COS               @                                 SIN               @                                 DBLE           
  @@                                                   
   @                                                   
   @                                                   D  @                              	                      D  @                              
                      
   @                                                   
   @                                                   
   @                                  
                
   @                                  
               
   @                                                 
    p          5 � p        r        5 � p        r                      #         @                                                     #EXIT_TPCORE%ALLOCATED                  @                                 ALLOCATED #         @                                                    #TPCORE_FVDAS%SUM    #DT    #AE    #IM    #JM    #KM    #JFIRST    #JLAST    #NG    #MG    #NQ    #AK    #BK    #U     #V !   #PS1 "   #PS2 #   #PS $   #Q %   #IORD &   #JORD '   #KORD (   #N_ADJ )   #XMASS *   #YMASS +   #FILL ,   #AREA_M2 -   #TCVV .   #WZ /                                      @                                 SUM           
   @                                  
                
   @                                  
                
  @@                                                  
  @@                                                  
  @@                                                  
   @                                                   
   @                                                   
   @                                                   
   @                                                   
   @                                                  
   @                                                 
    p           5 � p        r    n                                       1     5 � p        r    n                                      1                                    
   @                                                 
    p           5 � p        r    n                                       1     5 � p        r    n                                      1                                     
   @                                                 
              &                   &                   &                                                     
  @                             !                   
               &                   &                   &                                                    
D @@  �                          "                    
       5 � p        r    5 � p        r    p          5 � p        r      & 5 � p        r    5 � p        r        5 � p        r        5 � p        r    5 � p        r    p                                   
D @@  �                          #                    
       5 � p        r    5 � p        r    p          5 � p        r      & 5 � p        r    5 � p        r        5 � p        r        5 � p        r    5 � p        r    p                                   D  @  �                           $                    
       5 � p        r    5 � p        r    p          5 � p        r      & 5 � p        r    5 � p        r        5 � p        r        5 � p        r    5 � p        r    p                                    
D `@                            %                   
               &                   &                   &                                                     
  @@                              &                     
  @@                              '                     
  @@                              (                     
   @                              )                     
   @                            *                   
              &                   &                   &                                                     
   @                            +                   
              &                   &                   &                                                     
   @                              ,                    
  @@                            -                    
    p          5 � p        r        5 � p        r                                
   @                             .     
               D @@                             /                    
 %        p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      5 � p        r        5 � p        r      5 � p        r      5 � p        r                      #         @                                 0                  #FZPPM%SIGN 1   #FZPPM%ABS 2   #FZPPM%MIN 3   #FZPPM%MAX 4   #KLMT 5   #DELP1 6   #WZ =   #DQ1 B   #QQ1 C   #FZ D   #J1P E   #JU1_GL F   #J2_GL G   #ILO ;   #IHI :   #JULO 9   #JHI 8   #ILONG H   #IVERT I   #I1 A   #I2 @   #JU1 ?   #J2 >   #K1 7   #K2 <                 @                            1     SIGN               @                            2     ABS               @                            3     MIN               @                            4     MAX           
  @@                              5                    
   @  �                           6                    
 �       5 � p        r 7     5 � p        r 8   5 � p        r 9   p        5 � p        r 9     5 � p        r :   5 � p 
       r ;   p        5 � p 
       r ;     & 5 � p 
       r ;   5 � p        r :     & 5 � p        r 9   5 � p        r 8     & 5 � p        r 7   5 � p        r <         5 � p        r :   5 � p 
       r ;   p            5 � p        r 8   5 � p        r 9   p            5 � p        r <   5 � p        r 7   p                                   
   @  �                           =                    
 �       5 � p        r 7     5 � p        r >   5 � p        r ?   p        5 � p        r ?     5 � p        r @   5 � p        r A   p        5 � p        r A     & 5 � p        r A   5 � p        r @     & 5 � p        r ?   5 � p        r >     & 5 � p        r 7   5 � p        r <         5 � p        r @   5 � p        r A   p            5 � p        r >   5 � p        r ?   p            5 � p        r <   5 � p        r 7   p                                   
D  @  �                           B                    
 �        5 � p        r 7     5 � p        r 8   5 � p        r 9   p        5 � p        r 9     5 � p        r :   5 � p 
       r ;   p        5 � p 
       r ;     & 5 � p 
       r ;   5 � p        r :     & 5 � p        r 9   5 � p        r 8     & 5 � p        r 7   5 � p        r <         5 � p        r :   5 � p 
       r ;   p            5 � p        r 8   5 � p        r 9   p            5 � p        r <   5 � p        r 7   p                                   
   @  �                           C                    
 �       5 � p        r 7     5 � p        r 8   5 � p        r 9   p        5 � p        r 9     5 � p        r :   5 � p 
       r ;   p        5 � p 
       r ;     & 5 � p 
       r ;   5 � p        r :     & 5 � p        r 9   5 � p        r 8     & 5 � p        r 7   5 � p        r <         5 � p        r :   5 � p 
       r ;   p            5 � p        r 8   5 � p        r 9   p            5 � p        r <   5 � p        r 7   p                                   D  @  �                           D                    
 �        5 � p        r 7     5 � p        r 8   5 � p        r 9   p        5 � p        r 9     5 � p        r :   5 � p 
       r ;   p        5 � p 
       r ;     & 5 � p 
       r ;   5 � p        r :     & 5 � p        r 9   5 � p        r 8     & 5 � p        r 7   5 � p        r <         5 � p        r :   5 � p 
       r ;   p            5 � p        r 8   5 � p        r 9   p            5 � p        r <   5 � p        r 7   p                                    
   @                              E                     
   @                              F                     
   @                              G                     
   @                              ;                     
   @                              :                     
   @                              9                     
   @                              8                     
   @                              H                     
   @                              I                     
   @                              A                     
   @                              @                     
   @                              ?                     
   @                              >                     
   @                              7                     
   @                              <              �   .      fn#fn &   �   ;   b   uapp(TPCORE_FVDAS_MOD    	        INIT_TPCORE #     ?      INIT_TPCORE%REPEAT     U  <      INIT_TPCORE%COS     �  <      INIT_TPCORE%SIN !   �  =      INIT_TPCORE%DBLE    
  @   a   INIT_TPCORE%IM    J  @   a   INIT_TPCORE%JM    �  @   a   INIT_TPCORE%KM #   �  @   a   INIT_TPCORE%JFIRST "   
  @   a   INIT_TPCORE%JLAST    J  @   a   INIT_TPCORE%NG    �  @   a   INIT_TPCORE%MG    �  @   a   INIT_TPCORE%DT    
  @   a   INIT_TPCORE%AE !   J  �   a   INIT_TPCORE%CLAT    �  c       EXIT_TPCORE &   a  B      EXIT_TPCORE%ALLOCATED    �  q      TPCORE_FVDAS !     <      TPCORE_FVDAS%SUM     P  @   a   TPCORE_FVDAS%DT     �  @   a   TPCORE_FVDAS%AE     �  @   a   TPCORE_FVDAS%IM     	  @   a   TPCORE_FVDAS%JM     P	  @   a   TPCORE_FVDAS%KM $   �	  @   a   TPCORE_FVDAS%JFIRST #   �	  @   a   TPCORE_FVDAS%JLAST     
  @   a   TPCORE_FVDAS%NG     P
  @   a   TPCORE_FVDAS%MG     �
  @   a   TPCORE_FVDAS%NQ     �
  &  a   TPCORE_FVDAS%AK     �  &  a   TPCORE_FVDAS%BK      �   a   TPCORE_FVDAS%U    �  �   a   TPCORE_FVDAS%V !   �  �  a   TPCORE_FVDAS%PS1 !     �  a   TPCORE_FVDAS%PS2     �  �  a   TPCORE_FVDAS%PS       �   a   TPCORE_FVDAS%Q "   �  @   a   TPCORE_FVDAS%IORD "     @   a   TPCORE_FVDAS%JORD "   \  @   a   TPCORE_FVDAS%KORD #   �  @   a   TPCORE_FVDAS%N_ADJ #   �  �   a   TPCORE_FVDAS%XMASS #   �  �   a   TPCORE_FVDAS%YMASS "   T  @   a   TPCORE_FVDAS%FILL %   �  �   a   TPCORE_FVDAS%AREA_M2 "   H  @   a   TPCORE_FVDAS%TCVV     �  �  a   TPCORE_FVDAS%WZ      H      FZPPM    d  =      FZPPM%SIGN    �  <      FZPPM%ABS    �  <      FZPPM%MIN      <      FZPPM%MAX    U  @   a   FZPPM%KLMT    �    a   FZPPM%DELP1    �    a   FZPPM%WZ    �!    a   FZPPM%DQ1    �$    a   FZPPM%QQ1    �'    a   FZPPM%FZ    �*  @   a   FZPPM%J1P    9+  @   a   FZPPM%JU1_GL    y+  @   a   FZPPM%J2_GL    �+  @   a   FZPPM%ILO    �+  @   a   FZPPM%IHI    9,  @   a   FZPPM%JULO    y,  @   a   FZPPM%JHI    �,  @   a   FZPPM%ILONG    �,  @   a   FZPPM%IVERT    9-  @   a   FZPPM%I1    y-  @   a   FZPPM%I2    �-  @   a   FZPPM%JU1    �-  @   a   FZPPM%J2    9.  @   a   FZPPM%K1    y.  @   a   FZPPM%K2 