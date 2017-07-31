function rule_write ( order, filename, x, w, r )

%*****************************************************************************80
%
%% RULE_WRITE writes a quadrature rule to a file.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ORDER, the order of the rule.
%
%    Input, string FILENAME, specifies the output files.
%    write files 'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining 
%    weights, abscissas, and region.
%
%    Input, real X(ORDER), the abscissas.
%
%    Input, real W(ORDER), the weights.
%
%    Input, real R(2), the region.
%
  filename_x = strcat ( filename, '_x.txt' );
  filename_w = strcat ( filename, '_w.txt' );
  filename_r = strcat ( filename, '_r.txt' );

  fprintf ( 1, '\n' );
  fprintf ( 1,'  Creating quadrature files.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  "Root" file name is   "%s".\n', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Weight file will be   "%s".\n', filename_w );
  fprintf ( 1, '  Abscissa file will be "%s".\n', filename_x );
  fprintf ( 1, '  Region file will be   "%s".\n', filename_r );

  r8mat_write ( filename_w, 1, order, w' );
  r8mat_write ( filename_x, 1, order, x' );
  r8mat_write ( filename_r, 1, 2,     r' );

  return
end
