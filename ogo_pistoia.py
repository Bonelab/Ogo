"""
pistoia.py

Generate tables of standard post-processing quantities.

Copyright (c) 2010-2016, Numerics88 Solutions Ltd.
http://www.numerics88.com/
See LICENSE for details.
"""


import sys
#from N88ReportedError import N88ReportedError
from n88tools.N88ReportedError import N88ReportedError
import numpy
from numpy.core import *


def pistoia():

    import os
    import argparse
    from collections import OrderedDict
    import warnings
    import numbers
    from netCDF4 import Dataset

    global table_count

    # ------------------------------------------------------------------------
    # Parse options

    parser = argparse.ArgumentParser (
        description="""
    Calculate Pistoia yield critera.

    This is a method if estimating yield strength from linear solutions.
    In general, it would be preferable to use a non-linear elastoplastic
    model to calculate yield strengths. This utility is provided mostly
    for comparing with older results.
    """,
        )

    parser.add_argument ("--output_file", "-o",
        help="Output file. If not specified, output will go to STDOUT.")

    parser.add_argument ("--constraint", "-n", default="top_displacement",
        help="""Specify the constraint (i.e. boundary condition or applied load)
    to use for analysis. This is the surface to which forces are applied. The default
    ("top_displacement") will work for models generated with n88modelgenerator.""")

    parser.add_argument ("--rotation_center", "-c",
        type = lambda x: numpy.fromstring(x, dtype=float, sep=","),
        help="""Specify the spatial center used for calculation of angular
    quantities. The argument must be given as a triplet of coordinates.
    The default is read from the n88model.  If not specified, either on the
    command line or in the n88model file, no angular quantities will be
    calculated.""")

    parser.add_argument ("--crit_vol", "-cv", type = float, default = 7,
        help="Define the critical volume for the computed failure load.")

    parser.add_argument ("--crit_ees", "-cs", type = float, default = 0.011,
        help="Define the critical strain for the computed failure load.")

    parser_sets = parser.add_mutually_exclusive_group()

    parser_sets.add_argument ("--include", "-i",
        help="Only elements with the specified material IDs will be included in the calculation. Multiple IDs can be specified in a comma-delimited list (e.g. 100,101,105).")

    parser_sets.add_argument ("--exclude", "-e",
        help="All elements with the specified material IDs will be excluded from the calculation. Multiple IDs can be specified in a comma-delimited list (e.g. 100,101,105).")

    parser.add_argument ("model_file")

    args = parser.parse_args()
    crit_vol = args.crit_vol
    crit_ees = args.crit_ees

    if args.rotation_center != None:
        if len(args.rotation_center) != 3:
            raise N88ReportedError ("ERROR: rotation_center must be a triplet of values.")

    # For now this is hard-wired.
    args.twist_threshold = 1E-6

    if args.output_file == None:
        out = sys.stdout
    else:
        out = open (args.output_file, "wt")

    # ------------------------------------------------------------------------
    # Get information about the n88model file

    # Get the netCDF group handles from the n88model file.

    root = Dataset (args.model_file, "r")
    activeSolution = root.groups["Solutions"].groups[root.ActiveSolution]
    activeProblem = root.groups["Problems"].groups[activeSolution.Problem]
    activePart = root.groups["Parts"].groups[activeProblem.Part]
    hexahedrons = activePart.groups["Elements"].groups["Hexahedrons"]
    nodeValues = activeSolution.groups["NodeValues"]
    elementValues = activeSolution.groups["ElementValues"]
    materialTable = activePart.groups["MaterialTable"]
    materialDefinitions = root.groups["MaterialDefinitions"]
    nodeSetsGroup = root.groups["Sets"].groups["NodeSets"]
    elementSetsGroup = root.groups["Sets"].groups["ElementSets"]

    # Determine some constants to do with array sizes

    num = zeros (8,int)
    num[:] = hexahedrons.variables["NodeNumbers"][0] - 1
    coord_var = activePart.variables["NodeCoordinates"]
    spacing = array ([coord_var[num[1],0] - coord_var[num[0],0],
                     coord_var[num[2],1] - coord_var[num[0],1],
                     coord_var[num[4],2] - coord_var[num[0],2]])
    numberOfNodes = coord_var.shape[0]
    dimensionality = coord_var.shape[1]
    del num, coord_var
    numberOfElements = hexahedrons.variables["NodeNumbers"].shape[0]
    numberOfNodesPerElement = hexahedrons.variables["NodeNumbers"].shape[1]
    sizeOfMaterialTable = materialTable.variables["ID"].shape[0]
    nst = 6

    if args.rotation_center == None:
        try:
            args.rotation_center = activeProblem.RotationCenter
        except:
            args.rotation_center = None

    # Array data that it is convenient to read now and keep around.
    # Note that most array data will be read as needed and subsequently
    # discarded for memory efficiency.

    materialIds = materialTable.variables["ID"][:]
    elementMaterials = hexahedrons.variables["MaterialID"][:]

    constraints = root.groups["Constraints"]
    try:
        selected_constraint = constraints.groups[args.constraint]
    except:
        raise N88ReportedError ('ERROR: Specified constraint "%s" not found.' % args.constraint)
    set_indices = selected_constraint.variables["NodeNumber"][:]
    set_indices -= 1   # Convert to 0-indexed.


    # ------------------------------------------------------------------------
    # Some constants used for formatting

    width = 76
    table_delimiter = "="*width + "\n"
    section_delimiter = "-"*width + "\n"
    subtable_delimiter = "."*width + "\n"
    any_entry = "%%-32s%%%d" % (width-32)
    text_entry = any_entry + "s\n"
    integer_entry = any_entry + "d\n"


    # ------------------------------------------------------------------------
    # Output functions


    def CalculateStats (data, stats):
        if "median" in stats:
            # If sorting we make a copy so we don't modify the original data
            if isinstance(data, numpy.ma.MaskedArray):
                x = sort(data.compressed())
            else:
                x = sort(data)
        else:
            if isinstance(data, numpy.ma.MaskedArray):
                # If a MaskedArray, get a copy of the compressed array - may be faster
                x = data.compressed()
            else:
                x = data
        nan = float('nan')
        n = len(x)
        if n == 0:
            total = 0
            average = nan
            rms = nan
            std_dev = nan
            minimum = nan
            maximum = nan
            skewness = nan
            kurtosis = nan
            median = nan
            perc05 = nan
            perc25 = nan
            perc75 = nan
            perc95 = nan
        else:
            total = sum(x)
            average = total/n
            if "std_dev" in stats:
                std_dev = sqrt(sum((x-average)**2)/n)
            if "rms" in stats:
                rms = sqrt(sum(x**2)/n)
            if "minimum" in stats:
                minimum = min(x)
            if "maximum" in stats:
                maximum = max(x)
            # Might get warning if stddev is zero: ignore these
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if "skewness" in stats:
                    skewness = sum(((x-average)/std_dev)**3)/n
                if "kurtosis" in stats:
                    kurtosis = sum(((x-average)/std_dev)**4)/n - 3
            median = x[n//2]
            perc05 = x[n//20]
            perc25 = x[n//4]
            perc75 = x[(3*n)//4]
            perc95 = x[(19*n)//20]
        values = OrderedDict()
        for v in stats:
            values[v] = locals()[v]
        return values


    def BeginTable (name):
        global table_count
        table_count += 1
        out.write (
            "\n" + table_delimiter
          + "Table %d: %s\n" % (table_count, name)
          + section_delimiter)


    def StatsTable (data, stats, labels=None):
        if len(data.shape) == 1:
            num_col = 1
        else:
            num_col = data.shape[1]
        left_width = width - 11*6
        if labels:
            out.write (" "*left_width + ("%11s"*num_col) % tuple(labels) + "\n")
        int_row_format = "%%-%ds" % left_width + "%11d"*num_col + "\n"
        float_row_format = "%%-%ds" % left_width + "%11.3E"*num_col + "\n"
        if num_col == 1:
            values = CalculateStats(data, stats)
            for stat_name, stat_value in values.items():
                if isinstance(stat_value, numbers.Integral):
                    row_format = int_row_format
                else:
                    row_format = float_row_format
                out.write (row_format % (stat_name, stat_value))
        else:
            values = []
            for i in range(num_col):
                values.append (CalculateStats(data[:,i], stats))
            for v in list(values[0].keys()):
                row_values = []
                for i in range(num_col):
                    row_values.append (values[i][v])
                if isinstance(row_values[0], numbers.Integral):
                    row_format = int_row_format
                else:
                    row_format = float_row_format
                out.write (row_format % (tuple([v] + row_values)))
        return values


    # ------------------------------------------------------------------------
    # Generate a mask
    # Note that the mask is True (or value 1) for elements that are to be excluded.

    if args.exclude != None:
        ids = numpy.unique(numpy.fromstring(args.exclude, sep=",", dtype=int))
        mask = zeros((numberOfElements,), bool)
        for i in ids:
            mask += (elementMaterials == i)
    elif args.include != None:
        ids = numpy.unique(numpy.fromstring(args.include, sep=",", dtype=int))
        mask = ones((numberOfElements,), bool)
        for i in ids:
            mask *= (elementMaterials != i)
    else:
        mask = None


    # ------------------------------------------------------------------------
    # Calculations start here

    # Get average displacement of specified node set.
    displacements = zeros (nodeValues.variables["Displacement"].shape, float64)
    displacements[:] = nodeValues.variables["Displacement"][:]
    set_data = displacements[set_indices]
    u_ns1 = numpy.mean(set_data, axis=0)

    # Get total force on specified node set.
    forces = zeros (nodeValues.variables["ReactionForce"].shape, float64)
    forces[:] = nodeValues.variables["ReactionForce"][:]
    set_data = forces[set_indices]
    rf_ns1 = numpy.sum(set_data, axis=0)

    if args.rotation_center != None:

        # Calculate average angular rotation on specified node set.
        p = zeros (activePart.variables["NodeCoordinates"].shape, float64)
        p[:] = activePart.variables["NodeCoordinates"][:]
        p -= args.rotation_center
        p_dp = p + displacements
        twist = numpy.arctan2 (numpy.roll(p_dp,1,axis=1), numpy.roll(p_dp,2,axis=1)) \
              - numpy.arctan2 (numpy.roll(p,1,axis=1), numpy.roll(p,2,axis=1))
        # Correct for going over branch cut
        twist += 2*pi * (twist < pi)
        twist -= 2*pi * (twist > pi)
       # Mask out any points (initial or final) that are too close to the rotation center
       # in the rotation plane.
        mask = (numpy.roll(p_dp,1,axis=1)**2 + numpy.roll(p_dp,2,axis=1)**2
                   > args.twist_threshold**2)[:, numpy.newaxis]
        mask *= (numpy.roll(p,1,axis=1)**2 + numpy.roll(p,2,axis=1)**2
                   > args.twist_threshold**2)[:, numpy.newaxis]
        twist_masked = numpy.ma.MaskedArray (twist, numpy.logical_not(mask))
        set_data = twist_masked[set_indices]
        rot_ns1 = set_data.mean(axis=0)

        # Calculate average torque on speficifed node set.
        rf = zeros (nodeValues.variables["ReactionForce"].shape, float64)
        rf[:] = nodeValues.variables["ReactionForce"][:]
        torque = numpy.roll(rf,1,axis=1) * numpy.roll(p,2,axis=1) \
               - numpy.roll(rf,2,axis=1) * numpy.roll(p,1,axis=1)
        set_data = torque[set_indices]
        torque_ns1 = numpy.sum(set_data, axis=0)

    f = rf_ns1
    u = u_ns1
    if args.rotation_center != None:
        t = torque_ns1
        rot = rot_ns1

    fixed_critical_volume = crit_vol
    fixed_critical_ees = crit_ees

    # Get typical Modulus
    modulus = zeros(numberOfElements,float)
    matNames = materialTable.variables['MaterialName'][:]
    for id,matName in zip (materialIds, matNames):
        material = materialDefinitions.groups[matName]
        if hasattr(material, 'E'):
            if material.Type == "LinearIsotropic":
                modulus_value = material.E
            elif material.Type == "LinearOrthotropic":
                modulus_value = max(material.E)
            else:
                raise N88ReportedError ("ERROR: Unsupported material type for Pistoia calculation: %s" % material.Type)
            indices = numpy.nonzero(elementMaterials == id)
            modulus[indices] = modulus_value
        else:
            if material.Type != "LinearIsotropic":
                raise N88ReportedError ("ERROR: Unsupported material type for Pistoia calculation: %s" % material.Type)

            modulus_values = numpy.array(material.variables['E'][:])
            for i, E in enumerate(modulus_values):
                indices = numpy.nonzero(elementMaterials == i)
                modulus[indices] = E
    sed = zeros (elementValues.variables["StrainEnergyDensity"].shape, float64)
    sed[:] = elementValues.variables["StrainEnergyDensity"][:]
    if mask is not None:
        sed = numpy.ma.MaskedArray (sed, mask)
        modulus = numpy.ma.MaskedArray (modulus, mask)
    ees = sqrt(2*sed/modulus)
    del sed
    if mask is None:
        numberOfSelectedElements = numberOfElements
    else:
        numberOfSelectedElements = ees.count()
    if numberOfSelectedElements == 0:
        raise N88ReportedError ("ERROR: No elements in selection.")
    if mask is None:
        sort_indices = argsort (ees)
    else:
        sort_indices = ees.argsort(fill_value=inf)[:numberOfSelectedElements]
    ees = ees[sort_indices]
    el = int(numberOfSelectedElements*(100-fixed_critical_volume)/100)
    ees_at_crit_vol = ees[el]
    fixed_factor = fixed_critical_ees / ees_at_crit_vol
    failure_load = f * fixed_factor
    if args.rotation_center != None:
        torsional_failure_load = t * fixed_factor
    # Note that the use of numpy.float64 is to ensure IEEE math
    # (i.e. result=NAN on divide by zero instead of exception)
    # Might get warning if stddev is zero: ignore these
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stiffness = f/u
        if args.rotation_center != None:
            torsional_stiffness = t/rot

    table_count = 0
    BeginTable ("Pistoia Failure Load Estimate")
    out.write ("""Pistoia Criterion for Energy-Equivalent-Strain (EES)
    (*) Results valid for linear, isotropic model only.
    """)
    if args.rotation_center != None:
        out.write ("   (*) Warning: Torsion failure load not validated. Caution!\n")

    out.write (subtable_delimiter)
    float_entry = any_entry + ".3f\n"
    out.write (float_entry % ("Critical volume (%):", fixed_critical_volume))
    float_entry = "%%-40s%%%d" % (width-40) + ".4E\n"
    out.write (float_entry % ("Critical EES:", fixed_critical_ees))
    out.write (float_entry % ("EES at vol_crit:", ees_at_crit_vol))
    col_width = 12
    triple_entry = "%%-%ds" % (width-3*col_width) + ("%%%d.4E" % col_width)*3 + "\n"
    triple_entry_text = "%%-%ds" % (width-3*col_width) + ("%%%ds" % col_width)*3 + "\n"
    out.write (float_entry % ("Factor (from table):", fixed_factor))
    out.write (triple_entry % (("Failure load (RF * factor):",) + tuple(failure_load)))
    if args.rotation_center != None:
        out.write (triple_entry % (("Torsional failure load (T * factor):",) + tuple(torsional_failure_load)))
    out.write (section_delimiter)
    out.write ("Stiffness:\n")
    out.write (subtable_delimiter)
    out.write (triple_entry_text % ("", "x", "y", "z"))
    out.write (triple_entry % (("RF (node set 1):",) + tuple(f)))
    out.write (triple_entry % (("U (node set 1):",) + tuple(u)))
    out.write (triple_entry % (("Axial stiffness:",) + tuple(stiffness)))
    if args.rotation_center != None:
        out.write (triple_entry % (("T (node set 1):",) + tuple(t)))
        out.write (triple_entry % (("Rot (node set 1) [rad]:",) + tuple(rot)))
        out.write (triple_entry % (("Torsional stiffness:",) + tuple(torsional_stiffness)))

    out.write (section_delimiter)
    out.write ("Distribution of energy-equivalent-strain: EES = sqrt(2U/E).\n")
    out.write (subtable_delimiter)
    stats = ["average","std_dev","minimum","maximum","skewness","kurtosis","median"]
    StatsTable (ees, stats)

    out.write (section_delimiter)
    out.write ("Distribution of failed materials.\n")
    out.write (subtable_delimiter)
    g_mat_failed = elementMaterials[sort_indices[el:]]
    unique_materials = numpy.unique(elementMaterials)
    nMaterials = len(unique_materials)
    nels_failed = zeros (nMaterials, int)
    for i,m in enumerate(unique_materials):
        nels_failed[i] = sum(g_mat_failed == m)
    total = sum(nels_failed)
    #g_mat_failed = elementMaterials[sort_indices[el:]]
    #nels_failed = zeros (len(numpy.unique(elementMaterials)), int)
    #for i,m in enumerate(numpy.unique(elementMaterials)):
    #    nels_failed[i] = sum (g_mat_failed == m)
    #total = sum(nels_failed)
    col_width = 12
    left_margin = (width - 3*col_width)//2
    heading_format = " "*left_margin + ("%%%ds" % col_width)*3 + "\n"
    table_line = "-"*10
    row_format = " "*left_margin + "%%%ds" % col_width + "%%%dd" % col_width + "%%%d.2f" % col_width + "\n"
    out.write (heading_format % ("material","# els","%"))
    for i,m in enumerate(unique_materials):
        out.write (row_format % (m, nels_failed[i], numpy.float64(100)*nels_failed[i]/total))
    out.write (heading_format % ((table_line,)*3))
    out.write (row_format % ("total",total,100))

    out.write (section_delimiter)
    out.write ("Factor table:\n")
    crit_vol = arange(1,8,dtype=float)
    # crit_ees = array((0.005, 0.006, 0.007, 0.008, 0.009))
    crit_ees = array((0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.020))
    el = ((numberOfSelectedElements*(100-crit_vol)/100)).astype(int)
    factors = crit_ees[numpy.newaxis,:] / ees[el][:,numpy.newaxis]
    col_width = (width - 17)//5
    row_format = "%3d %10d |" + (" %%%d.3E" % col_width)*15 + "\n"
    heading_format = "%3s %10s |" + (" %%%ds" % col_width)*5 + "\n"
    out.write (heading_format % (("", "crit_vol", "crit_ees") + ("",)*4))
    heading_format = "%3s %10s |" + (" %%%ds" % col_width)*15 + "\n"
    out.write (heading_format % (("%","# vol") + tuple(crit_ees)))
    out.write (heading_format % (("--","-"*8) + ("-"*(col_width-2),)*15))
    for i in range(7):
        out.write (row_format % ((crit_vol[i], int(crit_vol[i]*numberOfSelectedElements))
                                 + tuple(factors[i])))

    out.write (table_delimiter)


def main():
    try:
        pistoia()
    except N88ReportedError as e:
        sys.stderr.write ("Error: " + e.message)
        sys.stderr.write ("\n")
        sys.exit (e.value)
    # Let other exceptions fall through to default python unhandled exception reporting.

if __name__ == "__main__":
    main()
