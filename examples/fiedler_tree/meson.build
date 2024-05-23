executable(
  'test_fiedler_tree',
  'test_fiedler_tree.c',
  dependencies : [butterfly_dep, m_dep]
)

foreach obj_file : [ 'dragon0.obj', 'dragon1.obj', 'dragon2.obj', 'dragon3.obj' ]
  fs.copyfile(obj_file, obj_file)
endforeach

scripts = ['plot_fiedler_node.py', 'plot_fiedler_perm.py']

foreach script : scripts
  fs.copyfile(script, script)
endforeach
