import yt
import numpy as np
import gc
gc.enable()

def most_significant_line(tree, tree_id, cosmo=None, z0_min=0.0):
    if cosmo == None:
        cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=tree.hubble_constant,
                                                 omega_lambda=tree.omega_lambda,
                                                 omega_matter=tree.omega_matter)

    # Find heads in the tree
    hnum = tree['uid'] == tree_id
    branch = tree[hnum][0]
    head_filter = (branch['tree', 'num_prog'] == 0) & (branch['tree', 'redshift'] >= z0_min)
    heads = branch['tree'][head_filter]

    # Set up t(z) lookup table
    sigfig = 3
    tz = {}

    # Search for the most "significant" line, meaning it has the most mass over time
    Mflux = tree.arr([0]*len(heads), 'Msun*s')
    for i, head in enumerate(heads):
        zround = round(head['redshift'], sigfig)
        if zround in tz:
            tlast = tz[zround]
        else:
            tlast = cosmo.t_from_z(head['redshift'])
            tz[zround] = tlast
        node = head.descendent
        while True:
            zround = round(node['redshift'], sigfig)
            if zround in tz:
                tnow = tz[zround]
            else:
                tnow = cosmo.t_from_z(node['redshift'])
                tz[zround] = tnow
            dt = tnow - tlast
            Mflux[i] += node['mass'] * dt
            #print (tz.keys(), zround)
            #print (head, node, node.descendent)
            #print (tnow, tlast, dt, Mflux[i], node['Mvir'], node['Mvir']*dt)
            tlast = tnow
            if hasattr(node, 'descendent'):
                node = node.descendent
                if node is None:
                    break
            else:
                break
    #for Mf, h in zip(Mflux, heads):
    #    print Mf, h['redshift'], h['Mvir']
    #print Mflux
    #print heads
    b_id = Mflux.argmax()
    newline = []
    node = heads[b_id]
    while True:
        newline.append(node)
        if hasattr(node, 'descendent'):
            node = node.descendent
            if node is None:
                break
        else:
            break

    # Create a dictionary of halo properties evolution arrays
    properties = {}
    N = len(newline)
    for f in tree.field_list + tree.derived_field_list:
        if isinstance(heads[b_id][f], yt.YTArray):
            unit = heads[b_id][f].units
        else:
            unit = 'dimensionless'
        Nd = heads[b_id][f].size
        if Nd > 1:
            properties[f] = tree.arr([[0]*Nd]*N, unit)
        else:
            properties[f] = tree.arr([0]*N, unit)
        for i, node in enumerate(newline):
            properties[f][i] = node[f]

    del Mflux
    del tz


    return newline, properties
