gm.plot               package:gmvalid               R Documentation

_P_l_o_t _g_r_a_p_h_i_c_a_l _m_o_d_e_l_s

_D_e_s_c_r_i_p_t_i_o_n:

     Plots given graphical models and writes provided measures next to
     the edges.

_U_s_a_g_e:

     gm.plot(model, significant = TRUE, data.analysis = 0)

_A_r_g_u_m_e_n_t_s:

   model: String vector with model formulas. See 'gm.modelsim'. 

significant: If TRUE only significant edges in the selected models are
          plotted (in solid lines). If FALSE also not significant edges
          are plotted as dashed lines. 

data.analysis: Upper-tri matrix with the measure for the edge between
          variables i and j (i>j) at matrix position [i,j]. If the
          length of the 'model' is bigger than 1, 'data.analysis' has
          to be a list of matrices. 

_D_e_t_a_i_l_s:

     The line width of the edges will depend on the size of the numbers
     in 'data.analysis'.

_V_a_l_u_e:

     TRUE

_N_o_t_e:

     Every use of the plot function opens a new window.

_A_u_t_h_o_r(_s):

     Fabian Sobotka, Marc Suling, Ronja Foraita 
      Bremen Institute for Prevention Research and Social Medicine 
      (BIPS)  <URL: http://www.bips.uni-bremen.de>

_S_e_e _A_l_s_o:

     'gm.analysis'

_E_x_a_m_p_l_e_s:

       gm.plot("ABC,CDE")

       gm.plot("VBA,EVC")
       
       gm.plot(c("ABC,CDE","AB,BC,CD,DE","ABC,DEF,GHI"))
       
       gm.plot("AB,AC",FALSE,matrix(0.5,nrow=3,ncol=3))

