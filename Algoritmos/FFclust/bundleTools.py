
#from soma import aims
import numpy as N
import os
import shutil
import math
#from soma import aims

def getBundleNames( bundlefile ):

  #get center names from bundle file
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  centers_num = len( bunlist ) / 2
  labels = []
  for i in range( centers_num ):

    ind = i * 2
    labels.append( bunlist[ ind ] )
  return labels

def getBundleNamesAndSizes( bundlefile ):

  #get center names from bundle file
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  curves_count = ns[ 'attributes' ][ 'curves_count' ]
  centers_num = len( bunlist ) / 2
  labels = []
  sizes = []
  for i in range( centers_num ):

    ind = i * 2
    labels.append( bunlist[ ind ] )
    prec_size = bunlist[ 1 ]
  for i in range( centers_num - 1 ):

    ind = i * 2 + 3
    prec_size_tmp = bunlist[ ind ]
    sizes.append( prec_size_tmp - prec_size )
    prec_size = prec_size_tmp
    sizes.append( curves_count - bunlist[ len( bunlist ) - 1 ] )
  return labels, sizes, curves_count

def getBundleSize( bundlefile ):

  #get center names from bundle file
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  size = ns[ 'attributes' ][ 'curves_count' ]
  return  size

def getBundleNb( bundlefile ):

  #get center names from bundle file
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  return  len( bunlist ) / 2

def allFibersToOneBundle( bundlefile, bundlename = '1', mode = 0 ):
  
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  curves_count = ns[ 'attributes' ][ 'curves_count' ]
  if ( mode == 1 ):

    bunlist = ns[ 'attributes' ][ 'bundles' ]
    bundlename = bunlist[ 0 ]
    bundles = '[ \'' + bundlename + '\', 0 ]'
    # write minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""
  open( bundlefile, 'w' ).write( minf % ( bundles, curves_count ) )

def changeBundleNameToNumber( bundlefile, bundleout, offset = 0 ):
  
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  curves_count = ns[ 'attributes' ][ 'curves_count' ]
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  # write minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""
  
  for i in range( len( bunlist ) / 2 ):
    ind = i * 2
    bunlist[ ind ] = str( i + offset )
    #print i, ' : len= ', len(points[i])

  bundlesstr = '['
  l = len( bunlist ) / 2
  for i in range( l - 1 ):
    ind = i * 2
    bundlesstr += ' \'' + bunlist[ ind ] + '\', ' \
                  + str( bunlist[ ind + 1 ] ) + ','

  bundlesstr += ' \'' + bunlist[ ind + 2 ] + '\', ' \
                + str( bunlist[ ind + 3 ] ) + ' ]'
  open( bundleout, 'w' ).write( minf % ( bundlesstr, curves_count ) )
  shutil.copyfile( bundlefile + 'data',  bundleout + 'data' )

def oneFiberPerBundle( bundlefile, bundleout, offset = 0 ):
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  curves_count = ns[ 'attributes' ][ 'curves_count' ]
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  # write minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""
  bunlist2 = []
  for i in range( curves_count ):

    bunlist2.append( str( i + offset ) )
    bunlist2.append( i )
    bundlesstr = '['
    l = len( bunlist2 ) / 2
  for i in range( l - 1 ):
    ind = i * 2
    bundlesstr += ' \'' + bunlist2[ ind ] + '\', '  \
                  + str(bunlist2[ ind + 1 ] ) + ','

  bundlesstr += ' \'' + bunlist2[ ind + 2 ] + '\', ' \
                + str( bunlist2[ ind + 3 ] ) + ' ]'
  open( bundleout, 'w' ).write( minf % ( bundlesstr, curves_count ) )
  shutil.copyfile( bundlefile + 'data',  bundleout + 'data')
56554516

def read_bundle( infile ):

  points = []
  bun_file = infile + 'data'
  os.path.getsize( bun_file )
  bytes = os.path.getsize( bun_file )
  num = bytes / 4

  num_count = 0
 # f = open( bun_file )
  with open(bun_file,'rb') as f:
    while num_count < num:
      # numero de puntos de la fibra

      p = N.frombuffer(f.read(4), 'i' )[ 0 ]
      vertex = N.frombuffer( f.read( p * 3 * 4 ), 'f' ).reshape( -1, 3 ) # lee coordenadas fibra
      points.append( vertex )
      #num_# commentunt = num_count + 1 + ( p * 3 )
      num_count = num_count + 1 + ( p * 3 )
#  f.close()

  #bundles = []
  return points#, bundles

def read_OneFiber( infile ):

  points = []
  bun_file = infile + 'data'
  os.path.getsize( bun_file )
  bytes = os.path.getsize( bun_file )
  num = bytes / 4
  f = open( bun_file )
  p = N.frombuffer( f.read( 4 ), 'i' )[ 0 ]
  vertex = N.frombuffer( f.read( p * 3 * 4 ), 'f' ).reshape( -1, 3 )
  points.append( vertex )
  f.close()

  #bundles = []
  return points#, bundles

def write_bundle( outfile, points ):

  #write bundles file
  f = open( outfile + 'data','wb' )
  ncount = len( points )
  for i in range( ncount ):
    f.write(N.array( [ len( points[ i ] ) ], N.int32 ).tostring() )
    f.write( points[ i ].ravel().tostring() )

  f.close()

  # wrtie minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""
  open( outfile, 'w' ).write(minf % ( [ 'points', 0 ], ncount ) )

def read_bundle_severalbundles( infile ):

  points = []
  bun_file = infile + 'data'
  os.path.getsize( bun_file )
  bytes = os.path.getsize( bun_file )
  num = bytes / 4

  ns = dict()
  exec(compile(open( infile ).read(), infile, 'exec'), ns)
  bundlescount = ns[ 'attributes' ][ 'bundles' ]
  curvescount = ns[ 'attributes' ][ 'curves_count' ]
  bunnames = []
  bunstart = []
  bun_num = len( bundlescount ) / 2
  count = 0
  for i in range( bun_num ):

    bunnames.append( bundlescount[ count ] )
    count = count + 1
    bunstart.append( bundlescount[ count ] )
    count = count + 1
    points.append( [] )

  bun_size = []
  if len( bunstart ) > 1:

    for b in range( len( bunstart ) - 1 ):

      bun_size.append( bunstart[ b + 1 ] - bunstart[ b ] )
      bun_size.append( curvescount - bunstart[ b + 1 ] )
  else:

    bun_size.append( curvescount )	
    num_count = 0
    f = open( bun_file )
    bun_count = 0
    num_count_bundle = 0
  while num_count < num:

    p = N.frombuffer( f.read( 4 ), 'i' )[ 0 ]
    vertex = N.frombuffer( f.read( p * 3 * 4 ), 'f' ).reshape( -1, 3 )
    points[ bun_count ].append( vertex )
    #print num_count, p, bun_count, num_count_bundle
    num_count_bundle = num_count_bundle + 1
    if num_count_bundle == bun_size[ bun_count ]:
      bun_count = bun_count + 1
      num_count_bundle = 0
      num_count = num_count + 1 + ( p * 3 )

  f.close()
  return points, bunnames


def write_bundle_severalbundles( outfile, points, bundles = [] ):

  #write bundles file
  f = open( outfile + 'data','wb' )
  ncount = 0
  for i in range( len( points ) ):

    size = len( points[ i ] )
    ncount = ncount + size
    bun = points[ i ]
    for i in range( size ):
      f.write( N.array( [ len( bun[ i ] ) ], N.int32 ).tostring() )
      f.write( bun[ i ].ravel().tostring() )
  f.close()
      # write minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""

  bundles_list = []
  ind = 0
  for i in range( len( points ) ):

    if bundles == []:

      bundles_list.append( str( i ) )
    else:

      bundles_list.append( bundles[ i ] )
    bundles_list.append( ind )
      #print( i, ' : len= ', len(points[i]))
    ind = ind + len( points[ i ] )

  bundlesstr = '['
  l = len( bundles_list ) // 2
  if l > 1:
    for i in range( l - 1 ):
      
      ind = i*2
      bundlesstr += ' \'' + bundles_list[ ind ] + '\', ' + str( bundles_list[ ind + 1 ] ) + ','
    bundlesstr += ' \'' + bundles_list[ ind+2 ] + '\', ' + str( bundles_list[ ind + 3 ] ) + ' ]'
  else:
    ind = 0
    bundlesstr += ' \'' + bundles_list[ ind ] + '\', ' + str( bundles_list[ ind + 1 ] ) + ']'
    
    

  open( outfile, 'w' ).write( minf % ( bundlesstr, ncount ) )


def getMinAndMaxFiberSteps(bunfile, returnDistances = False):
  points = read_bundle(bunfile)
  minv = 10000
  maxv = 0
  dists = []
  for p in points:
    for i in range(len(p)-1):
      p1 = p[i]
      p2 = p[i+1]
      x = p1[0]-p2[0]
      y = p1[1]-p2[1]
      z = p1[2]-p2[2]
      d = x*x + y*y + z*z
      if returnDistances:
        dists.append(math.sqrt(d))
      if d < minv:
        minv = d
      if d > maxv:
        maxv = d
        minv = math.sqrt(minv)
        maxv = math.sqrt(maxv)
  return minv, maxv, dists

def getSymmetricBundle(bunfile, bunout, onebundle_name = None):
  #if onebundle_name != None, all fibers are put into one bundle
  #named 'onebundle_name'
  points, bunnames = read_bundle_severalbundles(bunfile)
  points2 = []
  for bundle in points:
    bun = []
    for fiber in bundle:
      fiber2 = N.array(fiber)
      for p in fiber2:
        p[0] = -p[0]
        bun.append(fiber2)
        points2.append(bun)
        write_bundle_severalbundles(bunout, points2, bunnames)
  if onebundle_name != None:
    allFibersToOneBundle(bunout, onebundle_name)
  return

def bundleFlatHierarchy( bundlefile,
                         hiefile,
                         color_palette = 'colors256',
                         anat = None ):

  import anatomist.api as anatomist_api
  if anat == None:

    anat = anatomist_api.Anatomist()
    locals_var = locals()
    #read all bundles names
  bunlist = []
  ns = dict()
  exec(compile(open( bundlefile ).read(), bundlefile, 'exec'), ns)
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  bundles_num = len( bunlist ) / 2
  labels = []
  for i in range( bundles_num ):

    ind = i * 2
    labels.append( bunlist[ ind ] )

  #create hierarchy
  pl = anat.palettes()
  pal = pl.find( color_palette )
  color = pal.get()
  colorlist = []
  ind_color = 0
  h = aims.Hierarchy()
  h.setSyntax( 'hierarchy' )
  h['graph_syntax'] = 'RoiArg'
  all_tree = aims.Tree( True, 'fold_name' )
  all_tree[ 'name' ] = 'all'
  h.insert( all_tree )
  for i in range( len( labels ) ):

    indcol = ind_color % 256
    col = color.value( indcol )
    vcol = [ col.red(), col.green(), col.blue() ]
    tree = aims.Tree( True, 'fold_name' )
    tree[ 'color' ] = vcol
    tree[ 'name' ] = labels[ i ]
    all_tree.insert( tree )
    ind_color = ind_color + 1
    #store hierarchy
  aims.write( h, hiefile )

#creates a flat bundle hierarchy with different groups using a token of the bundle name
#example: if bundle_names = ['bundle_01_test1', 'bundle_02_test1', 'bundle_02_test2'], token_index = 1, delimiter = '_'
#          groups will be = ['01', '02' ]
def bundleFlatHierarchyByToken( bundlefile,
                                hiefile,
                                token_index = 0,
                                delimiter = '_',
                                color_palette = 'colors256',
                                anat = None ):

  import anatomist.api as anatomist_api
  if anat == None:

    anat = anatomist_api.Anatomist()

  #read all bundles names
  bunlist = []
  ns = dict()
  exec(compile(open(bundlefile).read(), bundlefile, 'exec'), ns)
  bunlist = ns[ 'attributes' ][ 'bundles' ]
  bundles_num = len( bunlist ) / 2
  labels = []
  for i in range( bundles_num ):

    ind = i * 2
    labels.append( bunlist[ ind ] )

  #create groups
  groups = dict()
  for lab in labels:

    group = lab.split( delimiter )[ token_index ]
    if ( group in groups) == False:

      groups[ group ] = []

    gr = groups[ group ]
    gr.append( lab )

  #create hierarchy
  pl = anat.palettes()
  pal = pl.find( color_palette )
  color = pal.get()
  colorlist = []
  ind_color = 0
  h = aims.Hierarchy()
  h.setSyntax( 'hierarchy' )
  h[ 'graph_syntax' ] = 'RoiArg'
  all_tree = aims.Tree( True, 'fold_name' )
  all_tree[ 'name' ] = 'all'
  h.insert( all_tree )
  for key in list(groups.keys()):

    tree = aims.Tree( True, 'fold_name' )
    tree[ 'name' ] = key
    all_tree.insert( tree )
    group = groups[ key ]

    for bun in group:

      indcol = ind_color % 256
      col = color.value( indcol)
      vcol = [ col.red(), col.green(), col.blue() ]
      ntree = aims.Tree( True, 'fold_name' )
      ntree[ 'color' ] = vcol
      ntree[ 'name' ] = bun
      tree.insert(ntree)
      ind_color = ind_color + 1

  #store hierarchy
  aims.write( h, hiefile )


#creates a flat bundle hierarchy with different groups using a list of bundle names
def bundleFlatHierarchyByNameList( names, labels,
                                   hiefile,
                                   color_palette = 'colors256',
                                   anat = None ):

  import anatomist.api as anatomist_api
  if anat == None:

    anat = anatomist_api.Anatomist()

  #create hierarchy
  pl = anat.palettes()
  pal = pl.find( color_palette )
  color = pal.get()
  colorlist = []
  ind_color = 0
  h = aims.Hierarchy()
  h.setSyntax( 'hierarchy' )
  h[ 'graph_syntax' ] = 'RoiArg'
  all_tree = aims.Tree( True, 'fold_name' )
  all_tree[ 'name' ] = 'all'
  h.insert( all_tree )
  labelnum = len(names)
  for i in range(labelnum):

    tree = aims.Tree( True, 'fold_name' )
    tree[ 'name' ] = names[i]
    all_tree.insert( tree )
    llist = labels[i]
    indcol = ind_color % 256
    col = color.value( indcol)
    vcol = [ col.red(), col.green(), col.blue() ]
    tree[ 'color' ] = vcol

    for bun in llist:
      ntree = aims.Tree( True, 'fold_name' )
      ntree[ 'name' ] = bun
      tree.insert(ntree)

    ind_color = ind_color + 1

  #store hierarchy
  aims.write( h, hiefile )

def getFiberLength(bunfile):
  points = read_bundle(bunfile)
  allDistances = []
  for p in points:
    dist=0
    for i in range(len(p)-1):
      p1 = p[i]
      p2 = p[i+1]
      x = p1[0]-p2[0]
      y = p1[1]-p2[1]
      z = p1[2]-p2[2]
      d = x*x + y*y + z*z
      dist+=math.sqrt(d)
      allDistances.append(dist)
  return allDistances

def getFiberMedianCentroid(bunfile, distmat, outfile):
  points=read_bundle(bunfile)
  dm_data=aims.read(distmat).arraydata()
  sumdist=dm_data[0].sum(axis=1)
  min_index=N.where(sumdist[0]==sumdist[0].min())[0][0]
  points_out=[]
  points_out.append(points[min_index])
  write_bundle(outfile,points_out)

def getBundleLongestFibers(bunfile,percentaje, outfile):
  points = read_bundle(bunfile)
  labels = getBundleNames(bunfile)
  lengths = getFiberLength(bunfile)
  data = []
  for i in range(len(lengths)):
    data.append([labels[i],lengths[i],points[i]])
    data.sort(key=lambda lis:lis[1], reverse=True)
    bundles=[]
    pointsOut=[]
    rang = int(round(len(data)*percentaje/100.))
    if rang >= 1:
      for i in range(rang):
        bundles.append(data[i][0])
        pointsOut.append([data[i][2]])
    else:
      for i in range(len(points)):
        bundles.append(data[i][0])
        pointsOut.append([data[i][2]])
        write_bundle_severalbundles(outfile,pointsOut,bundles)
        
def read_hie(filepath):
  data=open(filepath,'r')
  data_read=data.read()
  data.close()
  lines=data_read.split('\n')
  lines=[_f for _f in lines if _f]
  clean_lines=[a for a in lines if a.split(' ')[0] == 'name']
  hie=[]
  for line in clean_lines[1:]:
    if line.split(' ')[1]=='':
      new_cluster=[]
      cluster_index=clean_lines.index(line)+1
      while (cluster_index<len(clean_lines) and clean_lines[cluster_index].split(' ')[1]!=''):
        new_cluster.append(clean_lines[cluster_index].split(' ')[1])
        cluster_index+=1
        hie.append(new_cluster)
  return hie

    
   
