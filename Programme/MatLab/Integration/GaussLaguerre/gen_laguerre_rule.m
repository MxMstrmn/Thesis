function gen_laguerre_rule ( order, alpha, a, b, filename )

%*****************************************************************************80
%
%% GEN_LAGUERRE_RULE generates a Gauss-Laguerre rule.
%
%  Discussion:
%
%    This program computes a standard or exponentially weighted 
%    generalized Gauss-Laguerre quadrature rule and writes it to a file.
%
%    The user specifies:
%    * the ORDER (number of points) in the rule;
%    * ALPHA, the exponent of |X|;
%    * A, the left endpoint of integration;
%    * B, the scale factor in the exponential;
%    * FILENAME, the root name of the output files.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 February 2010
%
%  Author:
%
%    John Burkardt
%
%   time_stamp ( );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, 'GEN_LAGUERRE_RULE\n' );
%   fprintf ( 1, '  MATLAB version\n' );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  Compute a generalized Gauss-Laguerre rule for approximating\n' );
%   fprintf ( 1, '    Integral ( a <= x < oo ) |x-a|^ALPHA exp(-B*(x-a)) f(x) dx\n' );
%   fprintf ( 1, '  of order ORDER.\n' );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  The user specifies ORDER, ALPHA, A, B, and FILENAME.\n' );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  ORDER is the number of points.\n' );
%   fprintf ( 1, '  ALPHA is the exponent of |X|.\n' );
%   fprintf ( 1, '  A is the left endpoint (typically 0).\n' );
%   fprintf ( 1, '  B is the exponential scale factor (typically 1).\n' );
%   fprintf ( 1, '  FILENAME is used to generate 3 files:\n' );
%   fprintf ( 1, '  * filename_w.txt - the weight file\n' );
%   fprintf ( 1, '  * filename_x.txt - the abscissa file.\n' );
%   fprintf ( 1, '  * filename_r.txt - the region file.\n' );
%
%  Initialize the parameters.
%
  beta = 0.0;
%
%  Get ORDER.
%
  if ( nargin < 1 )
    order = input ( '  Enter the rule order ORDER.' );
  elseif ( ischar ( order ) )
    order = str2num ( order );
  end
%
%  Get ALPHA.
%
  if ( nargin < 2 )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  ALPHA is the exponent of |X| in the weighting function.\n' );
    fprintf ( 1, '  ALPHA is a real number strictly greater than -1.\n' );
    fprintf ( 1, '\n' );
    alpha = input ( '  Enter the value of ALPHA.' );
  elseif ( ischar ( alpha ) )
    alpha = str2num ( alpha );
  end
%
%  Get A.
%
  if ( nargin < 3 )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  A is the left endpoint, typically 0.\n' );
    fprintf ( 1, '\n' );
    a = input ( '  Enter the value of A.' );
  elseif ( ischar ( a ) )
    a = str2num ( a );
  end
%
%  Get B.
%
  if (  nargin < 4 )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  B is the exponential scale factor, typically 1.\n' );
    fprintf ( 1, '\n' );
    b = input ( '  Enter the value of B.' );
  elseif ( ischar ( b ) )
    b = str2num ( b );
  end
%
%  Get FILENAME:
%
  if ( nargin < 5 )
    fprintf ( 1,  '\n' );
    fprintf ( 1,  '  FILENAME specifies the ''root name'' of the quadrature files).\n' );
    filename = input ( '  Enter FILENAME as a quoted string:' );
  end
%
%  Input summary.
%
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  ORDER = %d\n', order );
%   fprintf ( 1, '  ALPHA = %f\n', alpha );
%   fprintf ( 1, '  A = %f\n', a );
%   fprintf ( 1, '  B = %f\n', b );
%   fprintf ( 1, '  FILENAME = "%s".\n', filename );
%
%  Construct the rule.
%
  kind = 5;
  [ x, w ] = cgqf ( order, kind, alpha, beta, a, b );
%
%  Write the rule.
%
  r = zeros ( 2, 1 );
  r(1) = a;
  r(2) = r8_huge ( );
  rule_write ( order, filename, x, w, r );
%
%  Terminate.
%
%   fprintf ( 1, '\n' );
%   fprintf ( 1, 'GEN_LAGUERRE_RULE:\n' );
%   fprintf ( 1, '  Normal end of execution.\n' );
%   fprintf ( 1, '\n' );
%   time_stamp ( );

  return
end