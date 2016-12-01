import fpl_utils
import orca


def get_corrected_gbl(halide, cation, ion="Pb"):
    """
    Retrun the corrected UMBO for simulations using the GBL solvent.

    **Parameters**

        halide: *list, str or str*
            A list of strings specifying what halide combination was used.
            If only one string is passed, then it is assumed to be uniform.
        cation: *str*
            The cation used.

    **Return**

        UMBO: *float*
            Corrected UMBO. Returns the average UMBO between the two oxygens
            of GBL.
    """
    name = fpl_utils.reduce_to_name(ion, halide, cation) + "_gbl_orca_1"
    data = orca.read(name)

    # Now that we have the data, we need to get the MBO's of interest
    MBOs_with_O = [atoms for atoms in data.MBO
                   if "O" in [a.element for a in atoms[0]]]
    O_indices = [a.index for atoms in MBOs_with_O for a in atoms[0]
                 if a.element == "O"]

    FBO = [1 if O_indices.count(i) == 2 else 2 for i in O_indices]
    MBOS = [m[1] for m in MBOs_with_O]
    UMBO = [fbo - mbo for fbo, mbo in zip(FBO, MBOS)]
    return sum(UMBO) / len(UMBO)


def get_corrected_nitromethane(halide, cation, ion="Pb"):
    """
    Retrun the corrected UMBO for simulations using the nitromethane solvent.

    **Parameters**

        halide: *list, str or str*
            A list of strings specifying what halide combination was used.
            If only one string is passed, then it is assumed to be uniform.
        cation: *str*
            The cation used.

    **Return**

        UMBO: *float*
            Corrected UMBO
    """
    name = fpl_utils.reduce_to_name(ion, halide, cation)
    name += "_nitromethane_orca_1"
    data = orca.read(name)

    # Now that we have the data, we need to get the MBO's of interest
    MBOs_with_O = [atoms for atoms in data.MBO
                   if "O" in [a.element for a in atoms[0]]]

    return sum([1.5 - m[1] for m in MBOs_with_O]) / len(MBOs_with_O)


correct_solvent = {"GBL": get_corrected_gbl,
                   "gbl": get_corrected_gbl,
                   "NITROMETHANE": get_corrected_nitromethane,
                   "nitromethane": get_corrected_nitromethane}

if __name__ == '__main__':
    UMBO = get_corrected_gbl("Cl", "FA")
    print "GBL: " + str(UMBO)
    UMBO = get_corrected_nitromethane("Cl", "FA")
    print "Nitromethane: " + str(UMBO)
