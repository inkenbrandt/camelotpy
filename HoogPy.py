from numpy import *

def HoogTransform(t, gamma, bigT, N, FUN, meth, init, d, work):
    '''Numerrical Laplace inversion using the Hoog method.
       t     (float variable) independent variable at which inverse
             function is to be evaluated
       gamma (float variable) parameter of inversion integral, governing
             accuracy of inversion
       bigT  (float variable) parameter used to discretize inversion
             integral, thus governing discretization error
       N     (integer variable) must be >= 2 and should be even; number
             of terms to be used in sum for approximation, thus
             governing truncation error
       FUN   name of user-supplied function which evaluates the complex
             values of transform fbar(p) for complex values of Laplace
             transform variable p
       meth  (integer variable) determines inversion method: 1=epsilon
             method; 2=ordinary quotient difference algorithm;
             3=acceleratedquotient difference algorithm
       init  (integer variable) =1 indicates this is the first call
             of the algorithm with these values of gamma, bigT, N and
             FUN. !=1 means array is not recomputed
       d     (complex array dimensioned in calling function) contains
             continued fraction coefficient if init >1 and meth=2 or 3
       work  (complex array dimensioned in calling function) work space
             used for calculating successive diagonals of the table

       This is a translation of the FORTRAN subroutine LAPADC originally by
       J.H. Knight. See code for more information.
       Requires numpy.'''
    #This is a translation of the FORTRAN subroutine LAPADC by J.H. Knight
    #of the Division of Environmental Mechanics, CSIRO, Canberra, ACT2601,
    #Australia originally coded in August 1986. The function performs a
    #numerical LaPlace inversion using the Hoog method. Two papers are
    #referenced in the original FORTRAN code:
    #DeHoog, F.R., Knight, J.H. and Stokes, A.N. An improved method for
    #numerical inversion of LaPlace transforms. SIAM J. SCI. STA. COMPUT.,
    # 3, PP. 357-366, 1982.
    #DeHoog, F.R., Knight, J.H. and Stokes, A.N. Subroutine LAPACC for LaPlace
    #inversion by acceleration of convergence of complex sum. ACM TRANS.
    #MATH. SOFTWARE, (manuscript in preperation).
    #Translated by Karl W. Bandilla in July 2009.
    #Comments are taken from the FORTRAN code

    #Following is the header of the original FORTRAN code

#      SUBROUTINE LAPADC(T,GAMMA,BIGT,N,F,FUN,METH,INIT,D,WORK)
#
#-----------------------------------------------------------------------
#
#     AUTHOR
#
#         J.H. KNIGHT
#         DIVISION OF ENVIRONMENTAL MECHANICS
#         CSIRO, CANBERRA, ACT 2601, AUSTRALIA
#         AUGUST, 1986
#
#     REFERENCES
#
#         DE HOOG, F.R., KNIGHT, J.H. AND STOKES, A.N.
#         AN IMPROVED METHOD FOR NUMERICAL INVERSION OF LAPLACE
#         TRANSFORMS. SIAM J. SCI. STAT. COMPUT., 3, PP. 357-366, 1982.
#
#         DE HOOG, F.R., KNIGHT, J.H. AND STOKES, A.N.
#         SUBROUTINE LAPACC FOR LAPLACE TRANSFORM INVERSION BY
#         ACCELERATION OF CONVERGENCE OF COMPLEX SUM.
#         ACM TRANS. MATH. SOFTWARE, (MANUSCRIPT IN PREPARATION).
#
#     ABSTRACT
#
#         LAPADC COMPUTES AN APPROXIMATION TO THE INVERSE LAPLACE
#         TRANSFORM F(T) CORRESPONDING TO THE LAPLACE TRANSFORM
#         FBAR(P), USING A RATIONAL APPROXIMATION TO THE FOURIER
#         SERIES RESULTING WHEN THE INVERSION INTEGRAL IS
#         DISCRETIZED USING THE TRAPEZOIDAL RULE.
#         THE RATIONAL APPROXIMATION IS A DIAGONAL PADE APPROXIMATION
#         COMPUTED USING EITHER THE EPSILON ALGORITHM, IN WHICH THE
#         EPSILON TABLE MUST BE RECOMPUTED AT EACH VALUE OF T, OR
#         USING THE QUOTIENT-DIFFERENCE ALGORITHM, IN WHICH THE
#         QUOTIENT-DIFFERENCE ALGORITHM IS USED TO COMPUTE THE
#         COEFFICIENTS OF THE ASSOCIATED CONTINUED FRACTION, AND THE
#         CONTINUED FRACTION IS EVALUATED AT EACH VALUE OF T.
#         IN ADDITION, THE EVALUATION OF THE CONTINUED FRACTION IS
#         MADE MORE ACCURATE WITH AN IMPROVED ESTIMATE OF THE
#         REMAINDER.
#
#     METHOD
#
#         A REAL VALUED INVERSE FUNCTION F(T) IS THE REAL PART OF THE
#         INTEGRAL
#                                INFINITY
#
#         (ONE/PI)*EXP(GAMMA*T)* INTEGRAL FBAR(GAMMA+I*W)*EXP(I*W*T)*DW
#
#                                W=ZERO
#
#         WITH FBAR THE LAPLACE TRANSFORM, GAMMA A REAL PARAMETER, AND
#         I**2=-1.  W IS THE VARIABLE OF INTEGRATION.
#         IF THIS IS DISCRETIZED USING THE TRAPEZOIDAL RULE WITH STEP
#         SIZE PI/BIGT, A FOURIER SERIES RESULTS, AND F(T) IS THE SUM
#         OF THE DISCRETIZATION ERROR
#
#               INFINITY
#
#               SUM     EXP(-TWO*GAMMA*K*BIGT)*F(T+TWO*K*BIGT)
#
#               K=1
#
#         AND THE REAL PART OF THE FOURIER SERIES
#
#                                    INFINITY
#
#           (ONE/BIGT)*EXP(GAMMA*T)* SUM  A(K)*EXP(I*K*PI*T/BIGT)
#
#                                    K=0
#
#         WITH THE COEFFICIENTS (A(K),K=0,INFINITY) GIVEN BY
#
#              A(0) = HALF*FBAR(GAMMA),
#
#              A(K) = FBAR(GAMMA+I*K*PI/BIGT),   K=1,INFINITY.
#
#         THIS EXPRESSION CAN ALSO BE CONSIDERED AS A POWER SERIES IN
#         THE VARIABLE Z=EXP(I*PI*T/BIGT), AND WE INTRODUCE A
#         TRUNCATION ERROR WHEN WE USE A FINITE POWER SERIES
#
#                        M2
#
#              V(Z,M2) = SUM  A(K)*Z**K,       M2 = M*2.
#
#                        K=0
#
#         THE EPSILON ALGORITHM CALCULATES THE DIAGONAL PADE
#         APPROXIMATION TO THIS FINITE POWER SERIES
#
#                        M                M
#
#              V(Z,M2) = SUM  B(K)*Z**K / SUM  C(K)*Z**K
#
#                        K=0              K=0
#
#                                    + HIGHER ORDER TERMS,   C(0)=ONE.
#
#         THE CALCULATION MUST BE DONE FOR EACH NEW VALUE OF T.
#         THIS METHOD IS USED FOR METH=1.
#
#         THE QUOTIENT-DIFFERENCE (QD) ALGORITHM IS ANOTHER WAY
#         OF CALCULATING THIS DIAGONAL PADE APPROXIMATION, IN THE FORM
#         OF A CONTINUED FRACTION
#
#         V(Z,M2) = D(0)/(ONE+D(1)*Z/(ONE+...+D(M2-1)*Z/(ONE+R2M(Z)))),
#
#         WHERE R2M(Z) IS THE UNKNOWN REMAINDER.
#         THE QD ALGORITHM GIVES THE COEFFICIENTS (D(K),K=0,M2), AND
#         THE CONTINUED FRACTION IS THEN EVALUATED BY FORWARD RECURRENCE,
#         FOR EACH NEW VALUE OF T, WITH R2M(Z) TAKEN AS D(M2)*Z.
#         THE ANSWERS SHOULD BE IDENTICAL TO THOSE GIVEN BY THE
#         EPSILON ALGORITHM, APART FROM ROUNDOFF ERROR.
#         THIS METHOD IS USED FOR METH=2.
#
#         IN ADDITION, THE CONVERGENCE OF THE CONTINUED FRACTION CAN
#         BE ACCELERATED WITH AN IMPROVED ESTIMATE OF THE REMAINDER,
#         SATISFYING
#
#               R2M + ONE + (D(M2-1)-D(M2))*Z = D(M2)*Z/R2M
#
#         THIS METHOD IS USED FOR METH=3.
#
#     USAGE
#
#         THE USER MUST CHOOSE THE PARAMETERS BIGT AND GAMMA WHICH
#        DETERMINE THE DISCRETIZATION ERROR, AND N(=M2) WHICH
#         DETERMINES THE TRUNCATION ERROR AND ROUNDOFF ERROR, GIVEN
#         THE CHOICE OF BIGT AND GAMMA.  AS N IS INCREASED, THE
#         TRUNCATION ERROR DECREASES BUT THE ROUNDOFF ERROR INCREASES,
#         SINCE THE EPSILON AND QD ALGORITHMS ARE NUMERICALLY UNSTABLE.
#         IT IS DESIRABLE TO USE DOUBLE PRECISION FOR THE CALCULATIONS,
#         UNLESS THE COMPUTER KEEPS AT LEAST 12 SIGNIFICANT DECIMAL
#         DIGITS IN SINGLE PRECISION.
#         METH=1 USES THE EPSILON ALGORITHM, METH=2 USES THE QUOTIENT-
#         DIFFERENCE ALGORITHM, AND METH=3 USES THE QD ALGORITHM WITH
#         IMPROVED ESTIMATE OF THE REMAINDER ON EVALUATION.
#
#         AN ESTIMATE OF THE OPTIMAL VALUE OF N FOR GIVEN BIGT AND
#         GAMMA CAN BE GOT BY VARYING N AND CHOOSING A VALUE SLIGHTLY
#         LESS THAN THAT WHICH CAUSES THE RESULTS FOR METH=1 TO VARY
#         SIGNIFICANTLY FROM THOSE FOR METH=2.
#
#         CALCULATIONS USING THE QD ALGORITHM ARE MOST EFFICIENT
#         AND MOST ACCURATE IF LAPADC IS CALLED WITH METH=3, AND WITH
#         INIT=1 ON THE FIRST CALL WITH NEW VALUES OF BIGT, GAMMA, N
#         OR FUN, AND WITH INIT=2 ON SUBSEQUENT CALLS WITH NEW VALUES
#         OF T AND THE OTHER PARAMETERS UNCHANGED.
#
#
#     DESCRIPTION OF ARGUMENTS
#
#         T      -  REAL (DOUBLE PRECISION) VARIABLE.  ON INPUT, VALUE
#                   OF INDEPENDENT VARIABLE AT WHICH INVERSE FUNCTION
#                   F(T) IS TO BE APPROXIMATED.  UNCHANGED ON OUTPUT.
#         GAMMA  -  REAL (DOUBLE PRECISION) VARIABLE.  ON INPUT,
#                   CONTAINS PARAMETER OF INVERSION INTEGRAL, GOVERNING
#                   ACCURACY OF INVERSION.  UNCHANGED ON OUTPUT.
#         BIGT   -  REAL (DOUBLE PRECISION) VARIABLE.  ON INPUT,
#                   CONTAINS PARAMETER USED TO DISCRETIZE INVERSION
#                   INTEGRAL. GOVERNS DISCRETIZATION ERROR.
#                   UNCHANGED ON OUTPUT.
#         N      -  INTEGER VARIABLE, SHOULD BE EVEN .GE. 2.  ON INPUT,
#                   CONTAINS NUMBER OF TERMS TO BE USED IN SUM FOR
#                   APPROXIMATION.  DETERMINES TRUNCATION ERROR.
#                   UNCHANGED ON OUTPUT.
#         F      -  REAL (DOUBLE PRECISION) VARIABLE.  ON OUTPUT,
#                   CONTAINS VALUE OF INVERSE FUNCTION F(T), EVALUATED
#                   AT T.
#         FUN    -  NAME OF USER-SUPPLIED SUBROUTINE WHICH EVALUATES
#                   (DOUBLE) COMPLEX VALUES OF TRANSFORM FBAR(P)
#                   FOR (DOUBLE) COMPLEX VALUES OF LAPLACE
#                   TRANSFORM VARIABLE P.  CALLED IF METH .EQ. 1 .OR.
#                   INIT .EQ. 1.  MUST BE DECLARED EXTERNAL IN
#                   (SUB)PROGRAM CALLING LAPADC.
#         METH   -  INTEGER VARIABLE, DETERMINING METHOD TO BE USED TO
#                   ACCELERATE CONVERGENCE OF SUM.
#                   ON INPUT, METH .EQ. 1 INDICATES THAT THE EPSILON
#                   ALGORITHM IS TO BE USED.
#                   METH .EQ. 2 INDICATES THAT THE ORDINARY QUOTIENT
#                   DIFFERENCE ALGORITHM IS TO BE USED.
#                   METH .EQ. 3 INDICATES THAT THE QUOTIENT DIFFERENCE
#                   ALGORITHM IS TO BE USED, WITH FURTHER ACCELERATION
#                   OF THE CONVERGENCE OF THE CONTINUED FRACTION.
#                   UNCHANGED ON OUTPUT.
#         INIT   -  INTEGER VARIABLE. ON INPUT, USED ONLY IF Q-D
#                   ALGORITHM IS SPECIFIED (METH .EQ. 2 OR 3).
#                   INIT .EQ. 1 INDICATES THAT THIS IS THE FIRST CALL
#                   OF THE ALGORITHM WITH THESE VALUES OF GAMMA, BIGT,
#                   N, FUN.  IN THIS CASE, THE VALUES OF THE CONTINUED
#                   FRACTION COEFFICIENTS D(0:N) ARE RECALCULATED, AND
#                   THE CONTINUED FRACTION IS EVALUATED.
#                   INIT .NE. 1 INDICATES THAT THE VALUES OF GAMMA,
#                   BIGT, N, AND FUN ARE UNCHANGED SINCE THE LAST CALL.
#                   IN THIS CASE, THE INPUT VALUES OF THE CONTINUED
#                   FRACTION COEFFICIENTS D(0:N) ARE USED, AND THE
#                   CONTINUED FRACTION IS EVALUATED.
#                   UNCHANGED ON OUTPUT.
#         D      -  (DOUBLE) COMPLEX ARRAY, OF VARIABLE DIMENSION (0:N).
#                   MUST BE DIMENSIONED IN CALLING PROGRAM.
#                   ON INPUT, IF METH .EQ. 2 OR 3 .AND. INIT .NE. 1,
#                   MUST CONTAIN CONTINUED FRACTION COEFFICIENTS
#                   OUTPUT ON A PREVIOUS CALL.
#                   ON OUTPUT, IF METH .EQ. 2 OR 3, CONTAINS CONTINUED
#                   FRACTION COEFFICIENTS CALCULATED ON THIS CALL
#                   (INIT .EQ. 1) OR AS SUPPLIED ON INPUT (INIT .NE. 1).
#         WORK   -  (DOUBLE) COMPLEX ARRAY OF VARIABLE DIMENSION (0:N).
#                   MUST BE DIMENSIONED IN CALLING PROGRAM.
#                   WORK SPACE USED FOR CALCULATING SUCCESSIVE DIAGONALS
#                   OF THE TABLE IF METH .EQ. 1 .OR. INIT .EQ. 1.
#                   ON OUTPUT, CONTAINS LAST CALCULATED DIAGONAL OF
#                   TABLE.
#
#     SUBROUTINE CALLED
#
#           FUN(P,FBAR) - USER-SUPPLIED SUBROUTINE WHICH EVALUATES
#                   (DOUBLE) COMPLEX VALUES OF TRANSFORM FBAR(P)
#                   FOR (DOUBLE) COMPLEX VALUES OF LAPLACE
#                   TRANSFORM VARIABLE P.  CALLED IF METH .EQ. 1 .OR.
#                   INIT .EQ. 1.  MUST BE DECLARED EXTERNAL IN
#                   (SUB)PROGRAM WHICH CALLS SUBROUTINE LAPADC.


    if N < 2:
        print 'Order too low (N at least 2)'
        return -1.0

    zeror = 0.0
    zero = 0.0 + 0.0j
    half = 0.5 + 0.0j
    one = 1.0 + 0.0j
#    work = zeros(N+1,dtype=complex128)
#    d = zeros(N+1,dtype=complex128)
    #make order multiple of 2
    M2 = (N / 2) * 2
    factor = pi / bigT
 #   print 'factor', factor
    argR = t * factor
    Z = cos(argR) + sin(argR) * 1j
    if meth == 1 or init == 1:
        #new: get all A's first
        ##########################
        #need to think of way to do this
        ###########################
        #calculate table if epsilon algorithm is to be used, or
        #if this is the first call to quotient-difference algorithm
        #with these parameters and inverse transform
        A = FUN(gamma + zeror * 1j)
        AOld = half * A
        if AOld == 0.0:
            return 0.0
        argI = factor
        A = FUN(gamma + argI * 1j)
        #initialize table entries
        if meth == 1:
            ztoj = Z
            term = A * ztoj
            work[0] = AOld + term
            work[1] = one / term
        else:
            d[0] = AOld
            work[0] = zero
            work[1] = A / AOld
            d[1] = -1.0 * work[1]
            AOld = A
        #calculate successive diagonals of table
        for j in range(2, N + 1):
            #initialize calculation of diagonal
            old2 = work[0]
            old1 = work[1]
            argI += factor
            A = FUN(gamma + argI * 1j)
         #   print 'A', A, A - AOld
            if A == 0:
                return 0.0
            if meth == 1:
                #calculate next term and sum of power series
                ztoj = Z * ztoj
                term = A * ztoj
                work[0] += term
                work[1] = one / term
            else:
                work[0] = zero
                work[1] = A / AOld
                AOld = A
        #    print work[0:10]
            #calculate diagonal using Rhombus rules
            for i in range(2, j+1):
                old3 = old2
                old2 = old1
                old1 = work[i]
        #        print 'old1', old1
                if meth == 1:
                    #epsilon algorithm rule
                    work[i] = old3 + one / (work[i-1] - old2)
                else:
                    #quotient difference algorithm rules
                    if i % 2 == 0:
                        #difference form
       #                 print i, 'up', old3, work[i-1], old2
                        if old2 == old3:
      #                      print 'ouch', work[i-1]
                            work[i] = work[i-1]
                        elif old3 + work[i-1] == 0j:
     #                       print 'pow'
                            work[i] = old2
                        else:
                            work[i] = old3 - old2 + work[i-1] 
                    else:
                        #quotient form
    #                    print i, 'down', old3, work[i-1], old2
                        work[i] = old3 * work[i-1] / old2
   #                 print 'work[i-1] old2', work[i-1], old2, old3
                   # print 'work', work[i]
            #save continued fraction coefficients
            if meth != 1:
                d[j] = -1.0 * work[j]
   #             print 'd[j]', d[j]
 #   print 'W', work
  #  print 'd', d
    if meth == 1:
        #result of epsilon algorithm computation
        Result0 = work[N].real
    else:
        #evaluate continued fraction
        #initialize recurrence relations
        AOld1 = d[0]
        AOld2 = d[0]
        BOld1 = one + d[1] * Z
        BOld2 = one
        #use recurrence relations
        for j in range(2,N+1):
            if meth == 3 and j == N+1:
                #further acceleration on last iteration if required
                h2m = half * (one + (d[N] - d[N+1]) * Z)
                print (one + d[N+1] * Z / h2m**2)
                r2m = -1.0 * h2m * \
                      (one - sqrt((one + d[N+1] * Z / h2m**2)))
                A = AOld1 + r2m * AOld2
                B = BOld1 + r2m * BOld2
            else:
               # print 'd', j, d[j]
                A = AOld1 + d[j] * Z * AOld2
                #print 'A', A
                AOld2 = AOld1
                AOld1 = A
                B = BOld1 + d[j] * Z * BOld2
                #print 'B', B
                BOld2 = BOld1
                BOld1 = B
              #  print 'other', A, B
        ctemp = A /B
        #result of quotient difference algorithm evaluation
        Result0 = ctemp.real
    #return required approximate inverse
 #   print 'hhog', exp(gamma * t), Result0, bigT 
    return exp(gamma * t) * Result0 / bigT          

