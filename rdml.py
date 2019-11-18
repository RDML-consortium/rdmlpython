#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import sys
import os
import re
import datetime
import zipfile
import tempfile
import argparse
import math
import warnings
import numpy as np
from lxml import etree as ET

# Fixme: Delete
import datetime as dt

class RdmlError(Exception):
    """Basic exception for errors raised by the RDML-Python library"""
    def __init__(self, message):
        Exception.__init__(self, message)
    pass


class secondError(RdmlError):
    """Just to have, not used yet"""
    def __init__(self, message):
        RdmlError.__init__(self, message)
    pass


def _get_first_child(base, tag):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The first child lxml node element found or None.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return node
    return None


def _get_first_child_text(base, tag):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The text of first child node element found or an empty string.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return node.text
    return ""


def _get_first_child_bool(base, tag, triple=True):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        triple: If True, None is returned if not found, if False, False

    Returns:
        The a bool value of tag or if triple is True None.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return _string_to_bool(node.text, triple)
    if triple is False:
        return False
    else:
        return None


def _get_step_sort_nr(elem):
    """Get the number of the step eg. for sorting.

    Args:
        elem: The node element. (lxml node)

    Returns:
        The a int value of the step node nr.
    """

    if elem is None:
        raise RdmlError('A step element must be provided for sorting.')
    ret = _get_first_child_text(elem, "nr")
    if ret == "":
        raise RdmlError('A step element must have a \"nr\" element for sorting.')
    return int(ret)


def _sort_list_int(elem):
    """Get the first element of the array as int. for sorting.

    Args:
        elem: The 2d list

    Returns:
        The a int value of the first list element.
    """

    return int(elem[0])


def _sort_list_float(elem):
    """Get the first element of the array as float. for sorting.

    Args:
        elem: The 2d list

    Returns:
        The a float value of the first list element.
    """

    return float(elem[0])


def _sort_list_digital_PCR(elem):
    """Get the first column of the list as int. for sorting.

    Args:
        elem: The list

    Returns:
        The a int value of the first list element.
    """

    arr = elem.split("\t")
    return int(arr[0]), arr[4]


def _string_to_bool(value, triple=True):
    """Translates a string into bool value or None.

    Args:
        value: The string value to evaluate. (string)
        triple: If True, None is returned if not found, if False, False

    Returns:
        The a bool value of tag or if triple is True None.
    """

    if value is None or value == "":
        if triple is True:
            return None
        else:
            return False
    if type(value) is bool:
        return value
    if type(value) is int:
        if value != 0:
            return True
        else:
            return False
    if type(value) is str:
        if value.lower() in ['false', '0', 'f', '-', 'n', 'no']:
            return False
        else:
            return True
    return


def _value_to_booldic(value):
    """Translates a string, list or dic to a dictionary with true/false.

    Args:
        value: The string value to evaluate. (string)

    Returns:
        The a bool value of tag or if triple is True None.
    """

    ret = {}
    if type(value) is str:
        ret[value] = True
    if type(value) is list:
        for ele in value:
            ret[ele] = True
    if type(value) is dict:
        for key, val in value.items():
            ret[key] = _string_to_bool(val, triple=False)
    return ret


def _get_first_child_by_pos_or_id(base, tag, by_id, by_pos):
    """Get a child element of the base node with a given tag and position or id.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        by_id: The unique id to search for. (string)
        by_pos: The position of the element in the list (int)

    Returns:
        The child node element found or raise error.
    """

    if by_id is None and by_pos is None:
        raise RdmlError('Either an ' + tag + ' id or a position must be provided.')
    if by_id is not None and by_pos is not None:
        raise RdmlError('Only an ' + tag + ' id or a position can be provided.')
    allChildren = _get_all_children(base, tag)
    if by_id is not None:
        for node in allChildren:
            if node.get('id') == by_id:
                return node
        raise RdmlError('The ' + tag + ' id: ' + by_id + ' was not found in RDML file.')
    if by_pos is not None:
        if by_pos < 0 or by_pos > len(allChildren) - 1:
            raise RdmlError('The ' + tag + ' position ' + by_pos + ' is out of range.')
        return allChildren[by_pos]


def _add_first_child_to_dic(base, dic, opt, tag):
    """Adds the first child element with a given tag to a dictionary.

    Args:
        base: The base node element. (lxml node)
        dic: The dictionary to add the element to (dictionary)
        opt: If false and id is not found in base, the element is added with an empty string (Bool)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The dictionary with the added element.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            dic[tag] = node.text
            return dic
    if not opt:
        dic[tag] = ""
    return dic


def _get_all_children(base, tag):
    """Get a list of all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A list with all child node elements found or an empty list.
    """

    ret = []
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            ret.append(node)
    return ret


def _get_all_children_id(base, tag):
    """Get a list of ids of all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A list with all child id strings found or an empty list.
    """

    ret = []
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            ret.append(node.get('id'))
    return ret


def _get_number_of_children(base, tag):
    """Count all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A int number of the found child elements with the id.
    """

    counter = 0
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            counter += 1
    return counter


def _check_unique_id(base, tag, id):
    """Find all child elements with a given group and check if the id is already used.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        id: The unique id to search for. (string)

    Returns:
        False if the id is already used, True if not.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            if node.get('id') == id:
                return False
    return True


def _create_new_element(base, tag, id):
    """Create a new element with a given tag and id.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag. (string)
        id: The unique id of the new element. (string)

    Returns:
        False if the id is already used, True if not.
    """

    if id is None or id == "":
        raise RdmlError('An ' + tag + ' id must be provided.')
    if not _check_unique_id(base, tag, id):
        raise RdmlError('The ' + tag + ' id "' + id + '" must be unique.')

    return ET.Element(tag, id=id)


def _add_new_subelement(base, basetag, tag, text, opt):
    """Create a new element with a given tag and id.

    Args:
        base: The base node element. (lxml node)
        basetag: Child elements group tag. (string)
        tag: Child elements own tag, to be created. (string)
        text: The text content of the new element. (string)
        opt: If true, the element is optional (Bool)

    Returns:
        Nothing, the base lxml element is modified.
    """

    if opt is False:
        if text is None or text == "":
            raise RdmlError('An ' + basetag + ' ' + tag + ' must be provided.')
        ET.SubElement(base, tag).text = text
    else:
        if text is not None and text != "":
            ET.SubElement(base, tag).text = text


def _change_subelement(base, tag, xmlkeys, value, opt, vtype, id_as_element=False):
    """Change the value of the element with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        value: The text content of the new element.
        opt: If true, the element is optional (Bool)
        vtype: If true, the element is optional ("string", "int", "float")
        id_as_element: If true, handle tag "id" as element, else as attribute

    Returns:
        Nothing, the base lxml element is modified.
    """

    # Todo validate values with vtype
    goodVal = value
    if vtype == "bool":
        ev = _string_to_bool(value, triple=True)
        if ev is None or ev == "":
            goodVal = ""
        else:
            if ev:
                goodVal = "true"
            else:
                goodVal = "false"

    if opt is False:
        if goodVal is None or goodVal == "":
            raise RdmlError('A value for ' + tag + ' must be provided.')

    if tag == "id" and id_as_element is False:
        if base.get('id') != goodVal:
            par = base.getparent()
            groupTag = base.tag.replace("{http://www.rdml.org}", "")
            if not _check_unique_id(par, groupTag, goodVal):
                raise RdmlError('The ' + groupTag + ' id "' + goodVal + '" is not unique.')
            base.attrib['id'] = goodVal
        return

    # Check if the tag already excists
    elem = _get_first_child(base, tag)
    if elem is not None:
        if goodVal is None or goodVal == "":
            base.remove(elem)
        else:
            elem.text = goodVal
    else:
        if goodVal is not None and goodVal != "":
            new_node = ET.Element(tag)
            new_node.text = goodVal
            place = _get_tag_pos(base, tag, xmlkeys, 0)
            base.insert(place, new_node)


def _get_or_create_subelement(base, tag, xmlkeys):
    """Get element with a given tag, if not present, create it.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)

    Returns:
        The node element with the tag.
    """

    # Check if the tag already excists
    if _get_first_child(base, tag) is None:
        new_node = ET.Element(tag)
        place = _get_tag_pos(base, tag, xmlkeys, 0)
        base.insert(place, new_node)
    return _get_first_child(base, tag)


def _remove_irrelevant_subelement(base, tag):
    """If element with a given tag has no children, remove it.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)

    Returns:
        The node element with the tag.
    """

    # Check if the tag already excists
    elem = _get_first_child(base, tag)
    if elem is None:
        return
    if len(elem) == 0:
        base.remove(elem)


def _move_subelement(base, tag, id, xmlkeys, position):
    """Change the value of the element with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        id: The unique id of the new element. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        position: the new position of the element (int)

    Returns:
        Nothing, the base lxml element is modified.
    """

    pos = _get_tag_pos(base, tag, xmlkeys, position)
    ele = _get_first_child_by_pos_or_id(base, tag, id, None)
    base.insert(pos, ele)


def _move_subelement_pos(base, tag, oldpos, xmlkeys, position):
    """Change the value of the element with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        oldpos: The unique id of the new element. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        position: the new position of the element (int)

    Returns:
        Nothing, the base lxml element is modified.
    """

    pos = _get_tag_pos(base, tag, xmlkeys, position)
    ele = _get_first_child_by_pos_or_id(base, tag, None, oldpos)
    base.insert(pos, ele)


def _get_tag_pos(base, tag, xmlkeys, pos):
    """Returns a position were to add a subelement with the given tag inc. pos offset.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        pos: The position relative to the tag elements (int)

    Returns:
        The int number of were to add the element with the tag.
    """

    count = _get_number_of_children(base, tag)
    offset = pos
    if pos is None or pos < 0:
        offset = 0
        pos = 0
    if pos > count:
        offset = count
    return _get_first_tag_pos(base, tag, xmlkeys) + offset


def _get_first_tag_pos(base, tag, xmlkeys):
    """Returns a position were to add a subelement with the given tag.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)

    Returns:
        The int number of were to add the element with the tag.
    """

    listrest = xmlkeys[xmlkeys.index(tag):]
    counter = 0
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") in listrest:
            return counter
        counter += 1
    return counter


def _writeFileInRDML(rdmlName, fileName, data):
    """Writes a file in the RDML zip, even if it existed before.

    Args:
        rdmlName: The name of the RDML zip file
        fileName: The name of the file to write into the zip
        data: The data string of the file

    Returns:
        Nothing, modifies the RDML file.
    """

    needRewrite = False

    if os.path.isfile(rdmlName):
        with zipfile.ZipFile(rdmlName, 'r') as RDMLin:
            for item in RDMLin.infolist():
                if item.filename == fileName:
                    needRewrite = True

    if needRewrite:
        tempFolder, tempName = tempfile.mkstemp(dir=os.path.dirname(rdmlName))
        os.close(tempFolder)

        # copy everything except the filename
        with zipfile.ZipFile(rdmlName, 'r') as RDMLin:
            with zipfile.ZipFile(tempName, mode='w', compression=zipfile.ZIP_DEFLATED) as RDMLout:
                RDMLout.comment = RDMLin.comment
                for item in RDMLin.infolist():
                    if item.filename != fileName:
                        RDMLout.writestr(item, RDMLin.read(item.filename))
                if data != "":
                    RDMLout.writestr(fileName, data)

        os.remove(rdmlName)
        os.rename(tempName, rdmlName)
    else:
        with zipfile.ZipFile(rdmlName, mode='a', compression=zipfile.ZIP_DEFLATED) as RDMLout:
            RDMLout.writestr(fileName, data)


def _linearRegression(x, y, useVal):
    """A function which calculates the slope and/or the intercept.

    Args:
        x: The numpy array of the cycles
        y: The numpy array that contains the values in the WoL
        useVal: The logical numpy array

    Returns:
        An array with the slope and intercept.
    """

    x2 = x * x
    xy = x * y
    sumx = np.nansum(x, axis=1)
    sumy = np.nansum(y, axis=1)
    sumx2 = np.nansum(x2, axis=1)
    sumxy = np.nansum(xy, axis=1)
    n = np.nansum(useVal, axis=1)

    ssx = sumx2 - (sumx * sumx) / n
    sxy = sumxy - (sumx * sumy) / n

    slope = sxy / ssx
    intercept = (sumy / n) - slope * (sumx / n)
    return [slope, intercept]


def _linearRegression2(x, y, z, vecExcluded):
    """A function which calculates the slope and/or the intercept.

    Args:
        x: The numpy array of the cycles
        y: The numpy array that contains the values in the WoL
        z: The logical numpy array

    Returns:
        An array with the slope and intercept.
    """

    x2 = x * x
    xy = x * y
    sumx = np.nansum(x, axis=1)
    sumy = np.nansum(y, axis=1)
    sumx2 = np.nansum(x2, axis=1)
    sumxy = np.nansum(xy, axis=1)
    n = np.nansum(z, axis=1)

    with np.errstate(divide='ignore', invalid='ignore'):
        ssx = sumx2 - (sumx * sumx) / n
        sxy = sumxy - (sumx * sumy) / n

        slope = sxy / ssx
        intercept = (sumy / n) - slope * (sumx / n)

    for i in range(len(vecExcluded)):
        if vecExcluded[i] or n[i] == 1 or np.isnan(slope[i]):
            slope[i] = 0
            intercept[i] = 0

    return [slope, intercept]


def _findLRstop(fluor, z):
    """A function which finds the stop cycle.

    Args:
        fluor: The array with the fluorescence values
        z: The row to work on

    Returns:
        An int with the stop cycle.
    """

    # First and Second Derivative values calculation

    # Shifted matrix of the raw data
    secondLast = 3
    while secondLast <= fluor.shape[1] and (np.isnan(fluor[z, secondLast - 1]) or np.isnan(fluor[z, secondLast - 2]) or np.isnan(fluor[z, secondLast - 3])):
            secondLast += 1

    fluorShift = np.roll(fluor[z], 1, axis=0)  # Shift to right - real position is -0.5
    fluorShift[0] = np.nan

    # Subtraction of the shifted matrix to the raw data
    firstDerivative = fluor[z] - fluorShift
    FDMcycles = np.nanargmax(firstDerivative, axis=0) + 1

    # Shifted matrix of the firstDerivative
    firstDerivativeShift = np.roll(firstDerivative, -1, axis=0)  # Shift to left
    firstDerivativeShift[-1] = np.nan
    # Subtraction of the firstDerivative values to the shifted matrix
    secondDerivative = firstDerivativeShift - firstDerivative

    # Data may have nan FixME: nan should be ignoread as not present
    if FDMcycles + 2 <= fluor.shape[1]:
        if (not np.isnan(fluor[z, FDMcycles - 1]) and
            not np.isnan(fluor[z, FDMcycles]) and
            not np.isnan(fluor[z, FDMcycles + 1]) and
            fluor[z, FDMcycles + 1] > fluor[z, FDMcycles] > fluor[z, FDMcycles - 1]):
            FDMcycles = FDMcycles + 2
    else:
        FDMcycles = fluor.shape[1]

    maxLoopSD = 0.0
    fstop = fluor.shape[1]

    for j in range(secondLast, FDMcycles):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            tempLoopSD = np.mean(secondDerivative[j - 2: j + 1], axis=0)
        # The > 0.000000000001 is to avoid float differences to the pascal version
        if not np.isnan(tempLoopSD) and (tempLoopSD - maxLoopSD) > 0.000000000001:
            maxLoopSD = tempLoopSD
            fstop = j
    if fstop + 2 >= fluor.shape[1]:
        fstop = fluor.shape[1]

    return fstop


def _findLRstart(fluor, z, fstop):
    """A function which finds the start cycle.

    Args:
        fluor: The array with the fluorescence values
        z: The row to work on
        fstop: The stop cycle

    Returns:
        An int with the stop cycle.
    """

    fstart = fstop - 1

    # Find the first value that is not NaN
    firstNotNaN = 1  # Cycles so +1 to array
    while np.isnan(fluor[z, firstNotNaN - 1]) and firstNotNaN < fstop:
        firstNotNaN += 1

    # fstart might be NaN, so shift it to the first value
    while fstart > firstNotNaN and np.isnan(fluor[z, fstart - 1]):
        fstart -= 1

    # As long as there are no NaN and new values are increasing
    while (fstart >= firstNotNaN and
           not np.isnan(fluor[z, fstart - 2]) and
           fluor[z, fstart - 2] <= fluor[z, fstart - 1]):
        fstart -= 1

    fstartfix = fstart
    if (not np.isnan(fluor[z, fstart]) and
        not np.isnan(fluor[z, fstart - 1]) and
        not np.isnan(fluor[z, fstop - 1]) and
        not np.isnan(fluor[z, fstop - 2])):
        startstep = np.log10(fluor[z, fstart]) - np.log10(fluor[z, fstart - 1])
        stopstep = np.log10(fluor[z, fstop - 1]) - np.log10(fluor[z, fstop - 2])
        if startstep > 1.1 * stopstep:
            fstartfix += 1

    return [fstart, fstartfix]


def _ttestslopes(fluor, samp, fstop, fstart):
    """A function which finds the start cycle.

    Args:
        fluor: The array with the fluorescence values
        samp: The row to work on
        fstop: The stop cycle
        fstart: The start cycle

    Returns:
        An array with [slopelow, slopehigh].
    """

    # Both start with full range
    loopStart = [fstart[samp], fstop[samp]]
    loopStop = [fstart[samp], fstop[samp]]

    # Now find the center ignoring nan
    while (loopStart[1] - loopStop[0]) > 1:
        loopStart[1] -= 1
        loopStop[0] += 1
        while (loopStart[1] - loopStop[0]) > 1 and np.isnan(fluor[samp, loopStart[1] - 1]):
            loopStart[1] -= 1
        while (loopStart[1] - loopStop[0]) > 1 and np.isnan(fluor[samp, loopStop[1] - 1]):
            loopStop[0] += 1

    # basic regression per group
    ssx = [0, 0]
    sxy = [0, 0]
    cslope = [0, 0]
    for j in range(0, 2):
        sumx = 0.0
        sumy = 0.0
        sumx2 = 0.0
        sumxy = 0.0
        nincl = 0.0
        for i in range(loopStart[j], loopStop[j] + 1):
            if not np.isnan(fluor[samp, i - 1]):
                sumx += i
                sumy += np.log10(fluor[samp, i - 1])
                sumx2 += i * i
                sumxy += i * np.log10(fluor[samp, i - 1])
                nincl += 1
        ssx[j] = sumx2 - sumx * sumx / nincl
        sxy[j] = sumxy - sumx * sumy / nincl
        cslope[j] = sxy[j] / ssx[j]

    return [cslope[0], cslope[1]]


def _GetMeanMax(fluor, ampgroup, vecSkipSample, vecNoPlateau):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        fluor: The array with the fluorescence values
        ampgroup: The target group id
        vecSkipSample: Skip the sample
        vecNoPlateau: Sample has no plateau

    Returns:
        An float with the max mean.
    """

    maxFlour = np.nanmax(fluor[:, -11:], axis=1)

    maxFlour[vecSkipSample] = np.nan
    maxFlour[vecNoPlateau] = np.nan

    # Ignore all nan slices, to fix them below
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        maxMean = np.nanmean(maxFlour)
    if maxMean == np.nan:
        maxMean = np.nanmax(maxFlour)

    return maxMean


def _GetMeanEff(tarGroup, vecTarget, pcreff, skipsample, noplateau, shortloglin, IsUsedInWoL):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        ampgroup: The target group id
        pcreff: The array with the PCR efficiencies
        skipsample: Skip the sample
        noplateau: Sample ha no plateau
        shortloglin: Sample has short loglin phase
        IsUsedInWoL: Sample is used in WoL - recalculated for samples in ampgroup

    Returns:
        An array with [getmeaneff, effvar, IsUsedInWoL].
    """

    cnt = 0
    sumeff = 0.0
    sumeff2 = 0.0

    for j in range(0, len(pcreff)):
        if tarGroup is None or tarGroup == vecTarget[j]:
            if (not (skipsample[j] or noplateau[j] or shortloglin[j])) and pcreff[j] > 1.0:
                IsUsedInWoL[j] = True
                cnt += 1
                sumeff += pcreff[j]
                sumeff2 += pcreff[j] * pcreff[j]
            else:
                IsUsedInWoL[j] = False

    if cnt > 1:
        getmeaneff = sumeff / cnt
        effvar = (sumeff2 - (sumeff * sumeff) / cnt) / (cnt - 1)
    else:
        getmeaneff = 1.0
        effvar = 100

    return [getmeaneff, effvar, IsUsedInWoL]


def _do_cascade(fluor, samp, uplim, lowlim, ftval):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        ampgroup: The target group id
        pcreff: The array with the PCR efficiencies
        skipsample: Skip the sample
        noplateau: Sample ha no plateau
        shortloglin: Sample has short loglin phase
        IsUsedInWoL: Sample is used in WoL

    Returns:
        An array with [getmeaneff, effvar, IsUsedInWoL].
    """

    pstart = 0
    pstop = 0
    fstop = _findLRstop(fluor, samp)
    [start, fstart] = _findLRstart(fluor, samp, fstop)

    stop = np.nanargmax(fluor[samp, fstart - 1:]) + fstart

    if fluor[samp, start - 1] > uplim or fluor[samp, stop - 1] < lowlim:
        notinwindow = True
        if fluor[samp, start - 1] > uplim:
            pstart = start
            pstop = start
        if fluor[samp, stop - 1] < lowlim:
            pstart = stop
            pstop = stop
    else:
        notinwindow = False
        # look for pstop
        if fluor[samp, stop - 1] < uplim:
            pstop = stop
        else:
            for i in range(stop, start, -1):
                if fluor[samp, i - 1] > uplim > fluor[samp, i - 2]:
                    pstop = i - 1
        # look for pstart
        if fluor[samp, fstart - 1] > lowlim:
            pstart = fstart
        else:
            for i in range(stop - 1, start - 1, -1):
                if fluor[samp, i] > lowlim > fluor[samp, i - 1]:
                    pstart = i + 1
   # print('find ' + str(samp) + ' pstart[samp]: ' + str(fluor[samp,pstart-1]) + ' pstop[samp]: ' + str(fluor[samp,pstop-1]))

    sumx = 0.0
    sumy = 0.0
    sumx2 = 0.0
    sumy2 = 0.0
    sumxy = 0.0
    nincl = 0.0
    ssx = 0.0
    ssy = 0.0
    sxy = 0.0
    for i in range(pstart, pstop + 1):
        if not np.isnan(fluor[samp, i - 1]):
            sumx += i
            sumy += np.log10(fluor[samp, i - 1])
            sumx2 += i * i
            sumy2 += np.log10(fluor[samp, i - 1]) * np.log10(fluor[samp, i - 1])
            sumxy += i * np.log10(fluor[samp, i - 1])
            nincl += 1

    if nincl > 0:  # Fixme: Should be 1???
        ssx = sumx2 - sumx * sumx / nincl
        ssy = sumy2 - sumy * sumy / nincl
        sxy = sumxy - sumx * sumy / nincl
    if ssx > 0.0 and nincl > 0.0:
        correl = sxy / np.sqrt(ssx * ssy)
        cslope = sxy / ssx
        cinterc = sumy / nincl - cslope * sumx / nincl
    else:
        correl = 0.0
        cslope = 0.0
        cinterc = 0.0

    if notinwindow:
        ninclu = 0
    else:
        ninclu = pstop - pstart + 1
    if cslope > 0:
        ctvals = (ftval - cinterc) / cslope
    else:
        ctvals = 0.0
    pcreff = np.power(10, cslope)
    nnulls = np.power(10, cinterc)

    return ctvals, pcreff, nnulls, ninclu, correl


def _do_all_cascade(fluor, tarGroup, vecTarget, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin, ftval, skipsample, noplateau):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        ampgroup: The target group id
        pcreff: The array with the PCR efficiencies
        skipsample: Skip the sample
        noplateau: Sample ha no plateau
        shortloglin: Sample has short loglin phase
        IsUsedInWoL: Sample is used in WoL

    Returns:
        An array with [getmeaneff, effvar, IsUsedInWoL].
    """

    if tarGroup is None:
        uplim = np.power(10, upwin[0])  # Fixme: No logs
        lowlim = np.power(10, lowwin[0])
    else:
        uplim = np.power(10, upwin[tarGroup])  # Fixme: No logs
        lowlim = np.power(10, lowwin[tarGroup])

    for z in range(0, fluor.shape[0]):
        if (tarGroup is None or tarGroup == vecTarget[z]) and not (skipsample[z] or noplateau[z]):
            ctvals[z], pcreff[z], nnulls[z], ninclu[z], correl[z] = _do_cascade(fluor, z, uplim, lowlim, ftval)

    return ctvals, pcreff, nnulls, ninclu, correl


def _GetMeanFluStop(fluor, tarGroup, vecTarget, fstop, skipsample, noplateau):
    # realStop = fstop[:, np.newaxis] - 1
    # vmeanmax = fluor[:, realStop[:, 0]]
    meanmax = 0.0
    maxfluor = 0.0000001
    cnt = 0
    if tarGroup is None:
        for j in range(0, fluor.shape[0]):
            if not skipsample[j]:
                if not noplateau[j]:
                    cnt += 1
                    meanmax += fluor[j, fstop[j] - 1]
                else:
                    for i in range(0, fluor.shape[1]):
                        if fluor[j, i] > maxfluor:
                            maxfluor = fluor[j, i]
    else:
        for j in range(0, fluor.shape[0]):
            if tarGroup == vecTarget[j] and not skipsample[j]:
                if not noplateau[j]:
                    cnt += 1
                    meanmax += fluor[j, fstop[j] - 1]
                else:
                    for i in range(0, fluor.shape[1]):
                        if fluor[j, i] > maxfluor:
                            maxfluor = fluor[j, i]

    if cnt > 0:
        meanmax = meanmax / cnt
    else:
        meanmax = maxfluor
    return meanmax


def _GetMaxFluStart(fluor, tarGroup, vecTarget, fstop, fstart, skipsample, noplateau):
    maxstart = -10.0
    if tarGroup is None:
        for j in range(0, fluor.shape[0]):
            if not skipsample[j]:
                if fluor[j, fstart[j] - 1] > maxstart:
                    maxstart = fluor[j, fstart[j]]
    else:
        for j in range(0, fluor.shape[0]):
            if tarGroup == vecTarget[j] and not skipsample[j]:
                if fluor[j, fstart[j] - 1] > maxstart:
                    maxstart = fluor[j, fstart[j]]

    return 0.999 * maxstart


def _ApplyLogWindow(tarGroup, uplim, foldwidth, upwin, lowwin, yaxismax, yaxismin):
    # Fixme: No truncation needed
    uplim = np.trunc(1000 * uplim) / 1000
    lowlim = np.trunc(1000 * (uplim - foldwidth)) / 1000

    uplim = np.minimum(uplim, yaxismax)
    uplim = np.maximum(uplim, yaxismin)
    lowlim = np.minimum(lowlim, yaxismax)
    lowlim = np.maximum(lowlim, yaxismin)

    # Fixme: No fix needed
    if lowlim < 0.0:
        lowlim += 0.001

    if tarGroup is None:
        upwin[0] = uplim
        lowwin[0] = lowlim
    else:
        upwin[tarGroup] = uplim
        lowwin[tarGroup] = lowlim

    return upwin, lowwin


def _GetLogStepStop(fluor, tarGroup, vecTarget, fstop, skipsample, noplateau):
    cnt = 0
    step = 0.0
    for j in range(0, fluor.shape[0]):
        if (tarGroup is None or tarGroup == vecTarget[j]) and not (skipsample[j] or noplateau[j]):
            cnt = cnt + 1
            step = step + (np.log10(fluor[j, fstop[j] - 1]) - np.log10(fluor[j, fstop[j] - 2]))
    if cnt > 0:
        step = step / cnt
    else:
        step = np.log10(1.8)
    return step


def _Set_WoL(fluor, tarGroup, vecTarget, PointsInWoL, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin, yaxismax, yaxismin, fstop, fstart, grpftval, vecSkipSample, vecNoPlateau, vecShortloglin, vecIsUsedInWoL):
    for i in range(0, fluor.shape[0]):
        if tarGroup == vecTarget[i]:
            if fstop[i] == fluor.shape[1]:
                vecNoPlateau[i] = True
            else:
                vecNoPlateau[i] = False

    skipgroup = False
    stepsize = 0.2  # was 0.5, smaller steps help in finding WoL
    vareffar = np.zeros(60, dtype=np.float64)
    uplimar = np.zeros(60, dtype=np.float64)
    foldwidthar = np.zeros(60, dtype=np.float64)

    maxlim = _GetMeanFluStop(fluor, tarGroup, vecTarget, fstop, vecSkipSample, vecNoPlateau)
    if maxlim > 0.0:
        maxlim = np.log10(maxlim)
    else:
        skipgroup = True
    minlim = _GetMaxFluStart(fluor, tarGroup, vecTarget, fstop, fstart, vecSkipSample, vecNoPlateau)
    if minlim > 0.0:
        minlim = np.log10(minlim)
    else:
        skipgroup = True

    thismeaneff = 1.0
    if not skipgroup:
        foldwidth = PointsInWoL * _GetLogStepStop(fluor, tarGroup, vecTarget, fstop, vecSkipSample, vecNoPlateau)
        upwin, lowwin = _ApplyLogWindow(tarGroup, maxlim, foldwidth, upwin, lowwin, yaxismax, yaxismin)
        print("Grp: " + str(tarGroup) + " upwin: " + str(upwin[tarGroup]) + " lowwin: " + str(lowwin[tarGroup]))
        tctvals, tpcreff, tnnulls, tninclu, tcorrel = _do_all_cascade(fluor, tarGroup, vecTarget, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin,
                                                                      grpftval, vecSkipSample, vecNoPlateau)
        thismeaneff, thisvareff, vecIsUsedInWoL = _GetMeanEff(tarGroup, vecTarget, tpcreff, vecSkipSample, vecNoPlateau, vecShortloglin,
                                                              vecIsUsedInWoL)
        print("thismeaneff: " + str(thismeaneff))
        if thismeaneff < 1.001:
            skipgroup = True

    if not skipgroup:
        foldwidth = np.log10(np.power(thismeaneff, PointsInWoL))
        k = -1
        maxvar = 0.0
        maxvarstep = -1
        lastuplim = 2 + maxlim
        while True:
            k += 1
            step = np.log10(thismeaneff)
            uplim = maxlim - k * stepsize * step
            if uplim < lastuplim:
                upwin, lowwin = _ApplyLogWindow(tarGroup, uplim, foldwidth, upwin, lowwin, yaxismax, yaxismin)
                tctvals, tpcreff, tnnulls, tninclu, tcorrel = _do_all_cascade(fluor, tarGroup, vecTarget, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin,
                                                                              grpftval, vecSkipSample, vecNoPlateau)
                thismeaneff, thisvareff, vecIsUsedInWoL = _GetMeanEff(tarGroup, vecTarget, tpcreff, vecSkipSample,
                                                                      vecNoPlateau, vecShortloglin,
                                                                      vecIsUsedInWoL)
                foldwidth = np.log10(np.power(thismeaneff, PointsInWoL))
                if foldwidth < 0.5:
                    foldwidth = 0.5  # to avoid width = 0 above fstop
                upwin, lowwin = _ApplyLogWindow(tarGroup, uplim, foldwidth, upwin, lowwin, yaxismax, yaxismin)
                tctvals, tpcreff, tnnulls, tninclu, tcorrel = _do_all_cascade(fluor, tarGroup, vecTarget, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin,
                                                                              grpftval, vecSkipSample, vecNoPlateau)
                thismeaneff, thisvareff, vecIsUsedInWoL = _GetMeanEff(tarGroup, vecTarget, tpcreff, vecSkipSample,
                                                                      vecNoPlateau, vecShortloglin,
                                                                      vecIsUsedInWoL)
                if thisvareff > 0.0:
                    vareffar[k] = np.sqrt(thisvareff) / thismeaneff
                else:
                    vareffar[k] = 0.0
                if thisvareff > maxvar:
                    maxvar = thisvareff
                    maxvarstep = k
                uplimar[k] = uplim
                foldwidthar[k] = foldwidth
                lastuplim = uplim
            else:
                thisvareff = 0.0
            if k >= 60 or uplim - foldwidth / (PointsInWoL / 2.0) < minlim or thisvareff < 0.00000000001:
                break

        # corrections: start
        if thisvareff < 0.00000000001:
            k -= 1  # remove window with vareff was 0.0
        i = -1
        while True:
            i += 1
            if vareffar[i] < 0.000001:
                break
        i -= 1  # i = number of valid steps

        minsmooth = vareffar[0]
        minstep = 0  # default top window

        # next 3 if conditions on i: added to correct smoothing
        if i == 0:
            minstep = 0

        if 0 < i < 4:
            n = -1
            while True:
                n += 1
                if vareffar[n] < minsmooth:
                    minsmooth = vareffar[n]
                    minstep = n
                if n == i:
                    break
        if i >= 4:
            n = 0
            while True:
                n += 1
                smoothvar = 0.0
                for m in range(n - 1, n + 2):
                    smoothvar = smoothvar + vareffar[m]
                smoothvar = smoothvar / 3.0
                if smoothvar < minsmooth:
                    minsmooth = smoothvar
                    minstep = n

                if n >= i - 1 or n > maxvarstep:
                    break
        # corrections: stop

        # Calculate the final values again
        upwin, lowwin = _ApplyLogWindow(tarGroup, uplimar[minstep], foldwidthar[minstep], upwin, lowwin, yaxismax, yaxismin)
        ctvals, pcreff, nnulls, ninclu, correl = _do_all_cascade(fluor, tarGroup, vecTarget, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin,
                                                                 grpftval, vecSkipSample, vecNoPlateau)
    return ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin, vecIsUsedInWoL


def _setWol(baselineCorrectedData, zeroBaselineCorrectedData, lowerLimit, upperLimit):
    # Function that sets the window of linearity : matWol is the matrix that contains only the log of
    # the values in the Wol (or 0 if the values are not in the WoL), and logicalMatrix contains 1 for
    # these values

    logicalMatrix2 = np.logical_and(np.less(np.full(baselineCorrectedData.shape, lowerLimit), zeroBaselineCorrectedData),
                                    np.less(zeroBaselineCorrectedData, np.full(baselineCorrectedData.shape, upperLimit)))

    vecCycle = np.arange(1, baselineCorrectedData.shape[1] + 1, dtype=np.int)

    valuesInWol = np.log10(baselineCorrectedData) * logicalMatrix2
    valuesInWol[np.isnan(valuesInWol)] = 0

    cyclesInWol = vecCycle * logicalMatrix2
    cyclesInWol[np.isnan(cyclesInWol)] = 0

    return [valuesInWol, logicalMatrix2, cyclesInWol]


def _numpyTwoAxisSave(var, fileName):
    with np.printoptions(precision=3, suppress=True):
        np.savetxt(fileName, var, fmt='%.12f', delimiter='\t', newline='\n')


class Rdml:
    """RDML-Python library
    
    The root element used to open, write, read and edit RDML files.
    
    Attributes:
        _rdmlData: The RDML XML object from lxml.
        _node: The root node of the RDML XML object.
    """

    def __init__(self, filename=None):
        """Inits an empty RDML instance with new() or load RDML file with load().

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """
        
        self._rdmlData = None
        self._rdmlFilename = None
        self._node = None
        if filename:
            self.load(filename)
        else:
            self.new()

    def __getitem__(self, key):
        """Returns data of the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """
        if key == "version":
            return self.version()
        if key in ["dateMade", "dateUpdated"]:
            return _get_first_child_text(self._node, key)
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["version", "dateMade", "dateUpdated"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["dateMade", "dateUpdated", "id", "experimenter", "documentation", "dye",
                "sample", "target", "thermalCyclingConditions", "experiment"]

    def new(self):
        """Creates an new empty RDML object with the current date.

        Args:
            self: The class self parameter.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        data = "<rdml version='1.2' xmlns:rdml='http://www.rdml.org' xmlns='http://www.rdml.org'>\n<dateMade>"
        data += datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data += "</dateMade>\n<dateUpdated>"
        data += datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data += "</dateUpdated>\n</rdml>"
        self.loadXMLString(data)
        return

    def load(self, filename):
        """Load an RDML file with decompression of rdml_data.xml or an XML file. Uses loadXMLString().

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """    

        if zipfile.is_zipfile(filename):
            self._rdmlFilename = filename
            zf = zipfile.ZipFile(filename, 'r')
            try:
                data = zf.read('rdml_data.xml').decode('utf-8')
            except KeyError:
                raise RdmlError('No rdml_data.xml in compressed RDML file found.')
            else:
                self.loadXMLString(data)
            finally:
                zf.close()
        else:
            with open(filename, 'r') as txtfile:
                data = txtfile.read()
                if data:
                    self.loadXMLString(data)
                else:
                    raise RdmlError('File format error, not a valid RDML or XML file.')

    def save(self, filename):
        """Save an RDML file with compression of rdml_data.xml.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to save to.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        elem = _get_or_create_subelement(self._node, "dateUpdated", self.xmlkeys())
        elem.text = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data = ET.tostring(self._rdmlData, pretty_print=True)
        _writeFileInRDML(filename, 'rdml_data.xml', data)

    def loadXMLString(self, data):
        """Create RDML object from xml string. !ENTITY and DOCSTRINGS will be removed.

        Args:
            self: The class self parameter.
            data: The xml string of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        # To avoid some xml attacs based on
        # <!ENTITY entityname "replacement text">
        data = re.sub(r"<\W*!ENTITY[^>]+>", "", data)
        data = re.sub(r"!ENTITY", "", data)
        try:
            self._rdmlData = ET.ElementTree(ET.fromstring(data.encode('utf-8')))
            # Change to bytecode and defused?
        except ET.XMLSyntaxError:
            raise RdmlError('XML load error, not a valid RDML or XML file.')
        self._node = self._rdmlData.getroot()
        if self._node.tag.replace("{http://www.rdml.org}", "") != 'rdml':
            raise RdmlError('Root element is not \'rdml\', not a valid RDML or XML file.')
        rdml_version = self._node.get('version')
        # Remainder: Update version in new() and validate()
        if rdml_version not in ['1.0', '1.1', '1.2', '1.3']:
            raise RdmlError('Unknown or unsupported RDML file version.')

    def validate(self, filename=None):
        """Validate the RDML object against its schema or load file and validate it.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            A string with the validation result as a two column table.
        """

        notes = ""
        if filename:
            try:
                vd = Rdml(filename)
            except RdmlError as err:
                notes += 'RDML file structure:\tFalse\t' + str(err) + '\n'
                return notes
            notes += "RDML file structure:\tTrue\tValid file structure.\n"
        else:
            vd = self
        version = vd.version()
        rdmlws = os.path.dirname(os.path.abspath(__file__))
        if version == '1.0':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        elif version == '1.3':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_3_CR.xsd'))
        else:
            notes += 'RDML version:\tFalse\tUnknown schema version' + version + '\n'
            return notes
        notes += "RDML version:\tTrue\t" + version + "\n"

        xmlschema = ET.XMLSchema(xmlschema_doc)
        result = xmlschema.validate(vd._rdmlData)
        if result:
            notes += 'Schema validation result:\tTrue\tRDML file is valid.\n'
        else:
            notes += 'Schema validation result:\tFalse\tRDML file is not valid.\n'
        log = xmlschema.error_log
        for err in log:
            notes += 'Schema validation error:\tFalse\t'
            notes += "Line %s, Column %s: %s \n" % (err.line, err.column, err.message)
        return notes

    def isvalid(self, filename=None):
        """Validate the RDML object against its schema or load file and validate it.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            True or false as the validation result.
        """

        if filename:
            try:
                vd = Rdml(filename)
            except RdmlError as err:
                return False
        else:
            vd = self
        version = vd.version()
        rdmlws = os.path.dirname(os.path.abspath(__file__))
        if version == '1.0':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        elif version == '1.3':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_3_CR.xsd'))
        else:
            return False
        xmlschema = ET.XMLSchema(xmlschema_doc)
        result = xmlschema.validate(vd._rdmlData)
        if result:
            return True
        else:
            return False

    def version(self):
        """Returns the version string of the RDML object.

        Args:
            self: The class self parameter.

        Returns:
            A string of the version like '1.1'.
        """

        return self._node.get('version')

    def migrate_version_1_0_to_1_1(self):
        """Migrates the rdml version from v1.0 to v1.1.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.0':
            raise RdmlError('RDML version for migration has to be v1.0.')

        exp = _get_all_children(self._node, "thirdPartyExtensions")
        if len(exp) > 0:
            ret.append("Migration to v1.1 deleted \"thirdPartyExtensions\" elements.")
        for node in exp:
            self._node.remove(node)

        hint = ""
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "react")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "data")
                    for node4 in exp4:
                        exp5 = _get_all_children(node4, "quantity")
                        for node5 in exp5:
                            hint = "Migration to v1.1 deleted react data \"quantity\" elements."
                            node4.remove(node5)
        if hint != "":
            ret.append(hint)

        xml_keys = ["description", "documentation", "xRef", "type", "interRunCalibrator",
                    "quantity", "calibratorSample", "cdnaSynthesisMethod",
                    "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]
        exp1 = _get_all_children(self._node, "sample")
        for node1 in exp1:
            hint = ""
            exp2 = _get_all_children(node1, "templateRNAQuantity")
            if len(exp2) > 0:
                templateRNAQuantity = _get_first_child_text(node1, "templateRNAQuantity")
                node1.remove(exp2[0])
                if templateRNAQuantity != "":
                    hint = "Migration to v1.1 modified sample \"templateRNAQuantity\" element without loss."
                    ele = _get_or_create_subelement(node1, "templateRNAQuantity", xml_keys)
                    _change_subelement(ele, "value", ["value", "unit"], templateRNAQuantity, True, "float")
                    _change_subelement(ele, "unit", ["value", "unit"], "ng", True, "float")
            if hint != "":
                ret.append(hint)
            hint = ""
            exp2 = _get_all_children(node1, "templateRNAQuantity")
            if len(exp2) > 0:
                templateDNAQuantity = _get_first_child_text(node1, "templateDNAQuantity")
                node1.remove(exp2[0])
                if templateDNAQuantity != "":
                    hint = "Migration to v1.1 modified sample \"templateDNAQuantity\" element without loss."
                    ele = _get_or_create_subelement(node1, "templateDNAQuantity", xml_keys)
                    _change_subelement(ele, "value", ["value", "unit"], templateDNAQuantity, True, "float")
                    _change_subelement(ele, "unit", ["value", "unit"], "ng", True, "float")
            if hint != "":
                ret.append(hint)

        xml_keys = ["description", "documentation", "xRef", "type", "amplificationEfficiencyMethod",
                    "amplificationEfficiency", "detectionLimit", "dyeId", "sequences", "commercialAssay"]
        exp1 = _get_all_children(self._node, "target")
        all_dyes = {}
        hint = ""
        for node1 in exp1:
            hint = ""
            dye_ele = _get_first_child_text(node1, "dyeId")
            node1.remove(_get_first_child(node1, "dyeId"))
            if dye_ele == "":
                dye_ele = "conversion_dye_missing"
                hint = "Migration to v1.1 created target nonsense \"dyeId\"."
            forId = _get_or_create_subelement(node1, "dyeId", xml_keys)
            forId.attrib['id'] = dye_ele
            all_dyes[dye_ele] = True
        if hint != "":
            ret.append(hint)
        for dkey in all_dyes.keys():
            if _check_unique_id(self._node, "dye", dkey):
                new_node = ET.Element("dye", id=dkey)
                place = _get_tag_pos(self._node, "dye", self.xmlkeys(), 999999)
                self._node.insert(place, new_node)

        xml_keys = ["description", "documentation", "experimenter", "instrument", "dataCollectionSoftware",
                    "backgroundDeterminationMethod", "cqDetectionMethod", "thermalCyclingConditions", "pcrFormat",
                    "runDate", "react"]
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                old_format = _get_first_child_text(node2, "pcrFormat")
                exp3 = _get_all_children(node2, "pcrFormat")
                for node3 in exp3:
                    node2.remove(node3)
                rows = "1"
                columns = "1"
                rowLabel = "ABC"
                columnLabel = "123"
                if old_format == "single-well":
                    rowLabel = "123"
                if old_format == "48-well plate; A1-F8":
                    rows = "6"
                    columns = "8"
                if old_format == "96-well plate; A1-H12":
                    rows = "8"
                    columns = "12"
                if old_format == "384-well plate; A1-P24":
                    rows = "16"
                    columns = "24"
                if old_format == "3072-well plate; A1a1-D12h8":
                    rows = "32"
                    columns = "96"
                    rowLabel = "A1a1"
                    columnLabel = "A1a1"
                if old_format == "32-well rotor; 1-32":
                    rows = "32"
                    rowLabel = "123"
                if old_format == "72-well rotor; 1-72":
                    rows = "72"
                    rowLabel = "123"
                if old_format == "100-well rotor; 1-100":
                    rows = "100"
                    rowLabel = "123"
                if old_format == "free format":
                    rows = "-1"
                    columns = "1"
                    rowLabel = "123"
                ele3 = _get_or_create_subelement(node2, "pcrFormat", xml_keys)
                _change_subelement(ele3, "rows", ["rows", "columns", "rowLabel", "columnLabel"], rows, True, "string")
                _change_subelement(ele3, "columns", ["rows", "columns", "rowLabel", "columnLabel"], columns, True, "string")
                _change_subelement(ele3, "rowLabel", ["rows", "columns", "rowLabel", "columnLabel"], rowLabel, True, "string")
                _change_subelement(ele3, "columnLabel", ["rows", "columns", "rowLabel", "columnLabel"], columnLabel, True, "string")
                if old_format == "48-well plate A1-F8" or \
                   old_format == "96-well plate; A1-H12" or \
                   old_format == "384-well plate; A1-P24":
                    exp3 = _get_all_children(node2, "react")
                    for node3 in exp3:
                        old_id = node3.get('id')
                        old_letter = ord(re.sub("\d", "", old_id).upper()) - ord("A")
                        old_nr = int(re.sub("\D", "", old_id))
                        newId = old_nr + old_letter * int(columns)
                        node3.attrib['id'] = str(newId)
                if old_format == "3072-well plate; A1a1-D12h8":
                    exp3 = _get_all_children(node2, "react")
                    for node3 in exp3:
                        old_id = node3.get('id')
                        old_left = re.sub("\D\d+$", "", old_id)
                        old_left_letter = ord(re.sub("\d", "", old_left).upper()) - ord("A")
                        old_left_nr = int(re.sub("\D", "", old_left)) - 1
                        old_right = re.sub("^\D\d+", "", old_id)
                        old_right_letter = ord(re.sub("\d", "", old_right).upper()) - ord("A")
                        old_right_nr = int(re.sub("\D", "", old_right))
                        newId = old_left_nr * 8 + old_right_nr + old_left_letter * 768 + old_right_letter * 96
                        node3.attrib['id'] = str(newId)
        self._node.attrib['version'] = "1.1"
        return ret

    def migrate_version_1_1_to_1_2(self):
        """Migrates the rdml version from v1.1 to v1.2.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.1':
            raise RdmlError('RDML version for migration has to be v1.1.')

        exp1 = _get_all_children(self._node, "sample")
        for node1 in exp1:
            hint = ""
            exp2 = _get_all_children(node1, "templateRNAQuality")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted sample \"templateRNAQuality\" element."
            if hint != "":
                ret.append(hint)
            hint = ""
            exp2 = _get_all_children(node1, "templateRNAQuantity")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted sample \"templateRNAQuantity\" element."
            if hint != "":
                ret.append(hint)
            hint = ""
            exp2 = _get_all_children(node1, "templateDNAQuality")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted sample \"templateDNAQuality\" element."
            if hint != "":
                ret.append(hint)
            hint = ""
            exp2 = _get_all_children(node1, "templateDNAQuantity")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted sample \"templateDNAQuantity\" element."
            if hint != "":
                ret.append(hint)
        self._node.attrib['version'] = "1.2"
        return ret

    def migrate_version_1_2_to_1_1(self):
        """Migrates the rdml version from v1.2 to v1.1.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.2':
            raise RdmlError('RDML version for migration has to be v1.2.')

        exp1 = _get_all_children(self._node, "sample")
        for node1 in exp1:
            hint = ""
            exp2 = _get_all_children(node1, "annotation")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.1 deleted sample \"annotation\" element."
            if hint != "":
                ret.append(hint)
            hint = ""
            exp2 = _get_all_children(node1, "templateQuantity")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.1 deleted sample \"templateQuantity\" element."
            if hint != "":
                ret.append(hint)

        exp1 = _get_all_children(self._node, "target")
        for node1 in exp1:
            hint = ""
            exp2 = _get_all_children(node1, "amplificationEfficiencySE")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.1 deleted target \"amplificationEfficiencySE\" element."
            if hint != "":
                ret.append(hint)

        hint = ""
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "react")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "data")
                    for node4 in exp4:
                        exp5 = _get_all_children(node4, "bgFluorSlp")
                        for node5 in exp5:
                            hint = "Migration to v1.1 deleted react data \"bgFluorSlp\" elements."
                            node4.remove(node5)
        if hint != "":
            ret.append(hint)

        self._node.attrib['version'] = "1.1"
        return ret

    def migrate_version_1_2_to_1_3(self):
        """Migrates the rdml version from v1.2 to v1.3.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.2':
            raise RdmlError('RDML version for migration has to be v1.2.')

        self._node.attrib['version'] = "1.3"
        return ret

    def migrate_version_1_3_to_1_2(self):
        """Migrates the rdml version from v1.3 to v1.2.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.3':
            raise RdmlError('RDML version for migration has to be v1.3.')

        hint = ""
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "react")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "partitions")
                    for node4 in exp4:
                        hint = "Migration to v1.2 deleted react \"partitions\" elements."
                        node3.remove(node4)
                    # No data element, no react element in v 1.2
                    exp5 = _get_all_children(node3, "data")
                    if len(exp5) == 0:
                        hint = "Migration to v1.2 deleted run \"react\" elements."
                        node2.remove(node3)

        if hint != "":
            ret.append(hint)

        self._node.attrib['version'] = "1.2"
        return ret

    def recreate_lost_ids(self):
        """Searches for lost ids and repairs them.

        Args:
            self: The class self parameter.

        Returns:
            A string with the modifications.
        """

        mess = ""

        # Find lost dyes
        foundIds = {}
        allTar = _get_all_children(self._node, "target")
        for node in allTar:
            forId = _get_first_child(node, "dyeId")
            if forId is not None:
                foundIds[forId.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "dye")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_dye(id=used_id, newposition=0)
                mess += "Recreated new dye: " + used_id + "\n"
        # Find lost thermalCycCon
        foundIds = {}
        allSam = _get_all_children(self._node, "sample")
        for node in allSam:
            subNode = _get_first_child(node, "cdnaSynthesisMethod")
            if subNode is not None:
                forId = _get_first_child(node, "thermalCyclingConditions")
                if forId is not None:
                    foundIds[forId.attrib['id']] = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                forId = _get_first_child(subNode, "thermalCyclingConditions")
                if forId is not None:
                    foundIds[forId.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "thermalCyclingConditions")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_therm_cyc_cons(id=used_id, newposition=0)
                mess += "Recreated thermal cycling conditions: " + used_id + "\n"
        # Find lost experimenter
        foundIds = {}
        allTh = _get_all_children(self._node, "thermalCyclingConditions")
        for node in allTh:
            subNodes = _get_all_children(node, "experimenter")
            for subNode in subNodes:
                foundIds[subNode.attrib['id']] = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                lastNodes = _get_all_children(subNode, "experimenter")
                for lastNode in lastNodes:
                    foundIds[lastNode.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "experimenter")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_experimenter(id=used_id, firstName="unknown first name", lastName="unknown last name", newposition=0)
                mess += "Recreated experimenter: " + used_id + "\n"
        # Find lost documentation
        foundIds = {}
        allSam = _get_all_children(self._node, "sample")
        for node in allSam:
            subNodes = _get_all_children(node, "documentation")
            for subNode in subNodes:
                foundIds[subNode.attrib['id']] = 0
        allTh = _get_all_children(self._node, "target")
        for node in allTh:
            subNodes = _get_all_children(node, "documentation")
            for subNode in subNodes:
                foundIds[subNode.attrib['id']] = 0
        allTh = _get_all_children(self._node, "thermalCyclingConditions")
        for node in allTh:
            subNodes = _get_all_children(node, "documentation")
            for subNode in subNodes:
                foundIds[subNode.attrib['id']] = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "documentation")
            for subNode in subNodes:
                foundIds[subNode.attrib['id']] = 0
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                lastNodes = _get_all_children(subNode, "documentation")
                for lastNode in lastNodes:
                    foundIds[lastNode.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "documentation")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_documentation(id=used_id, newposition=0)
                mess += "Recreated documentation: " + used_id + "\n"
        # Find lost sample
        foundIds = {}
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    lastNodes = _get_all_children(reactNode, "sample")
                    for lastNode in lastNodes:
                        foundIds[lastNode.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "sample")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_sample(id=used_id, type="unkn", newposition=0)
                mess += "Recreated sample: " + used_id + "\n"
        # Find lost target
        foundIds = {}
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "tar")
                        for lastNode in lastNodes:
                            foundIds[lastNode.attrib['id']] = 0
                    partNodes = _get_all_children(reactNode, "partitions")
                    for partNode in partNodes:
                        dataNodes = _get_all_children(partNode, "data")
                        for dataNode in dataNodes:
                            lastNodes = _get_all_children(dataNode, "tar")
                            for lastNode in lastNodes:
                                foundIds[lastNode.attrib['id']] = 0
        # Search in Table files
        if self._rdmlFilename is not None and self._rdmlFilename != "":
            if zipfile.is_zipfile(self._rdmlFilename):
                zf = zipfile.ZipFile(self._rdmlFilename, 'r')
                for item in zf.infolist():
                    if re.search("^partitions/", item.filename):
                        fileContent = zf.read(item.filename).decode('utf-8')
                        newlineFix = fileContent.replace("\r\n", "\n")
                        tabLines = newlineFix.split("\n")
                        header = tabLines[0].split("\t")
                        for cell in header:
                            if cell != "":
                                foundIds[cell] = 0
                zf.close()

        presentIds = []
        exp = _get_all_children(self._node, "target")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_target(id=used_id, type="toi", newposition=0)
                mess += "Recreated target: " + used_id + "\n"
        return mess

    def repair_rdml_file(self):
        """Searches for known errors and repairs them.

        Args:
            self: The class self parameter.

        Returns:
            A string with the modifications.
        """

        mess = ""
        mess += self.fixExclFalse()
        mess += self.fixDuplicateReact()

        return mess

    def fixExclFalse(self):
        """Searches in experiment-run-react-data for excl=false and deletes the elements.

        Args:
            self: The class self parameter.

        Returns:
            A string with the modifications.
        """

        mess = ""
        count = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "excl")
                        for lastNode in lastNodes:
                            if lastNode.text.lower() == "false":
                                count += 1
                                dataNode.remove(lastNode)

        if count > 0:
            mess = "The element excl=false was removed " + str(count) + " times!\n"

        return mess

    def fixDuplicateReact(self):
        """Searches in experiment-run-react for duplicates and keeps only the first.

        Args:
            self: The class self parameter.

        Returns:
            A string with the modifications.
        """

        mess = ""
        foundIds = {}
        count = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    tId = reactNode.attrib['id']
                    if tId not in foundIds:
                        foundIds[tId] = 0
                    else:
                        count += 1
                        subNode.remove(reactNode)

        if count > 0:
            mess = str(count) + " duplicate react elements were removed!\n"

        return mess

    def rdmlids(self):
        """Returns a list of all rdml id elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all rdml id elements.
        """

        exp = _get_all_children(self._node, "id")
        ret = []
        for node in exp:
            ret.append(Rdmlid(node))
        return ret

    def new_rdmlid(self, publisher, serialNumber, MD5Hash=None, newposition=None):
        """Creates a new rdml id element.

        Args:
            self: The class self parameter.
            publisher: Publisher who created the serialNumber (required)
            serialNumber: Serial Number for this file provided by publisher (required)
            MD5Hash: A MD5 hash for this file (optional)
            newposition: The new position of the element in the list (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = ET.Element("id")
        _add_new_subelement(new_node, "id", "publisher", publisher, False)
        _add_new_subelement(new_node, "id", "serialNumber", serialNumber, False)
        _add_new_subelement(new_node, "id", "MD5Hash", MD5Hash, True)
        place = _get_tag_pos(self._node, "id", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_rdmlid(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "id", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "id", None, oldposition)
        self._node.insert(pos, ele)

    def get_rdmlid(self, byposition=None):
        """Returns an experimenter element by position or id.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Rdmlid(_get_first_child_by_pos_or_id(self._node, "id", None, byposition))

    def delete_rdmlid(self, byposition=None):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "id", None, byposition)
        self._node.remove(elem)

    def experimenters(self):
        """Returns a list of all experimenter elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all experimenter elements.
        """

        exp = _get_all_children(self._node, "experimenter")
        ret = []
        for node in exp:
            ret.append(Experimenter(node))
        return ret

    def new_experimenter(self, id, firstName, lastName, email=None, labName=None, labAddress=None, newposition=None):
        """Creates a new experimenter element.

        Args:
            self: The class self parameter.
            id: Experimenter unique id
            firstName: Experimenters first name (required)
            lastName: Experimenters last name (required)
            email: Experimenters email (optional)
            labName: Experimenters lab name (optional)
            labAddress: Experimenters lab address (optional)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "experimenter", id)
        _add_new_subelement(new_node, "experimenter", "firstName", firstName, False)
        _add_new_subelement(new_node, "experimenter", "lastName", lastName, False)
        _add_new_subelement(new_node, "experimenter", "email", email, True)
        _add_new_subelement(new_node, "experimenter", "labName", labName, True)
        _add_new_subelement(new_node, "experimenter", "labAddress", labAddress, True)
        place = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_experimenter(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Experimenter unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "experimenter", id, self.xmlkeys(), newposition)

    def get_experimenter(self, byid=None, byposition=None):
        """Returns an experimenter element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Experimenter(_get_first_child_by_pos_or_id(self._node, "experimenter", byid, byposition))

    def delete_experimenter(self, byid=None, byposition=None):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "experimenter", byid, byposition)
        self._node.remove(elem)

    def documentations(self):
        """Returns a list of all documentation elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all documentation elements.
        """

        exp = _get_all_children(self._node, "documentation")
        ret = []
        for node in exp:
            ret.append(Documentation(node))
        return ret

    def new_documentation(self, id, text=None, newposition=None):
        """Creates a new documentation element.

        Args:
            self: The class self parameter.
            id: Documentation unique id
            text: Documentation descriptive test (optional)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "documentation", id)
        _add_new_subelement(new_node, "documentation", "text", text, True)
        place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_documentation(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Documentation unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "documentation", id, self.xmlkeys(), newposition)

    def get_documentation(self, byid=None, byposition=None):
        """Returns an documentation element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Documentation(_get_first_child_by_pos_or_id(self._node, "documentation", byid, byposition))

    def delete_documentation(self, byid=None, byposition=None):
        """Deletes an documentation element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "documentation", byid, byposition)
        self._node.remove(elem)

    def dyes(self):
        """Returns a list of all dye elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all dye elements.
        """

        exp = _get_all_children(self._node, "dye")
        ret = []
        for node in exp:
            ret.append(Dye(node))
        return ret

    def new_dye(self, id, description=None, newposition=None):
        """Creates a new dye element.

        Args:
            self: The class self parameter.
            id: Dye unique id
            description: Dye descriptive test (optional)
            newposition: Dye position in the list of dyes (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "dye", id)
        _add_new_subelement(new_node, "dye", "description", description, True)
        place = _get_tag_pos(self._node, "dye", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_dye(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Dye unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "dye", id, self.xmlkeys(), newposition)

    def get_dye(self, byid=None, byposition=None):
        """Returns an dye element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Dye(_get_first_child_by_pos_or_id(self._node, "dye", byid, byposition))

    def delete_dye(self, byid=None, byposition=None):
        """Deletes an dye element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "dye", byid, byposition)
        self._node.remove(elem)

    def samples(self):
        """Returns a list of all sample elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all sample elements.
        """

        exp = _get_all_children(self._node, "sample")
        ret = []
        for node in exp:
            ret.append(Sample(node))
        return ret

    def new_sample(self, id, type, newposition=None):
        """Creates a new sample element.

        Args:
            self: The class self parameter.
            id: Sample unique id (required)
            type: Sample type (required)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """

        if type not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
            raise RdmlError('Unknown or unsupported sample type value "' + type + '".')
        new_node = _create_new_element(self._node, "sample", id)
        _add_new_subelement(new_node, "sample", "type", type, False)
        place = _get_tag_pos(self._node, "sample", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_sample(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Sample unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "sample", id, self.xmlkeys(), newposition)

    def get_sample(self, byid=None, byposition=None):
        """Returns an sample element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Sample(_get_first_child_by_pos_or_id(self._node, "sample", byid, byposition))

    def delete_sample(self, byid=None, byposition=None):
        """Deletes an sample element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "sample", byid, byposition)
        self._node.remove(elem)

    def targets(self):
        """Returns a list of all target elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all target elements.
        """

        exp = _get_all_children(self._node, "target")
        ret = []
        for node in exp:
            ret.append(Target(node, self._rdmlFilename))
        return ret

    def new_target(self, id, type, newposition=None):
        """Creates a new target element.

        Args:
            self: The class self parameter.
            id: Target unique id (required)
            type: Target type (required)
            newposition: Targets position in the list of targets (optional)

        Returns:
            Nothing, changes self.
        """

        if type not in ["ref", "toi"]:
            raise RdmlError('Unknown or unsupported target type value "' + type + '".')
        new_node = _create_new_element(self._node, "target", id)
        _add_new_subelement(new_node, "target", "type", type, False)
        place = _get_tag_pos(self._node, "target", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_target(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Target unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "target", id, self.xmlkeys(), newposition)

    def get_target(self, byid=None, byposition=None):
        """Returns an target element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Target(_get_first_child_by_pos_or_id(self._node, "target", byid, byposition), self._rdmlFilename)

    def delete_target(self, byid=None, byposition=None):
        """Deletes an target element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "target", byid, byposition)
        self._node.remove(elem)

    def therm_cyc_cons(self):
        """Returns a list of all thermalCyclingConditions elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all target elements.
        """

        exp = _get_all_children(self._node, "thermalCyclingConditions")
        ret = []
        for node in exp:
            ret.append(Therm_cyc_cons(node))
        return ret

    def new_therm_cyc_cons(self, id, newposition=None):
        """Creates a new thermalCyclingConditions element.

        Args:
            self: The class self parameter.
            id: ThermalCyclingConditions unique id (required)
            newposition: ThermalCyclingConditions position in the list of ThermalCyclingConditions (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = _create_new_element(self._node, "thermalCyclingConditions", id)
        step = ET.SubElement(new_node, "step")
        ET.SubElement(step, "nr").text = "1"
        ET.SubElement(step, "lidOpen")
        place = _get_tag_pos(self._node, "thermalCyclingConditions", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_therm_cyc_cons(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: ThermalCyclingConditions unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "thermalCyclingConditions", id, self.xmlkeys(), newposition)

    def get_therm_cyc_cons(self, byid=None, byposition=None):
        """Returns an thermalCyclingConditions element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Therm_cyc_cons(_get_first_child_by_pos_or_id(self._node, "thermalCyclingConditions", byid, byposition))

    def delete_therm_cyc_cons(self, byid=None, byposition=None):
        """Deletes an thermalCyclingConditions element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "thermalCyclingConditions", byid, byposition)
        self._node.remove(elem)

    def experiments(self):
        """Returns a list of all experiment elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all experiment elements.
        """

        exp = _get_all_children(self._node, "experiment")
        ret = []
        for node in exp:
            ret.append(Experiment(node, self._rdmlFilename))
        return ret

    def new_experiment(self, id, newposition=None):
        """Creates a new experiment element.

        Args:
            self: The class self parameter.
            id: Experiment unique id (required)
            newposition: Experiment position in the list of experiments (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = _create_new_element(self._node, "experiment", id)
        place = _get_tag_pos(self._node, "experiment", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_experiment(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Experiments unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "experiment", id, self.xmlkeys(), newposition)

    def get_experiment(self, byid=None, byposition=None):
        """Returns an experiment element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Experiment(_get_first_child_by_pos_or_id(self._node, "experiment", byid, byposition), self._rdmlFilename)

    def delete_experiment(self, byid=None, byposition=None):
        """Deletes an experiment element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "experiment", byid, byposition)
        experiment = Experiment(elem, self._rdmlFilename)

        # Required to delete digital files
        runs = _get_all_children(elem, "run")
        for node in runs:
            run = Run(node, self._rdmlFilename)
            experiment.delete_run(byid=run["id"])

        # Now delete the experiment element
        self._node.remove(elem)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        allRdmlids = self.rdmlids()
        rdmlids = []
        for elem in allRdmlids:
            rdmlids.append(elem.tojson())

        allExperimenters = self.experimenters()
        experimenters = []
        for exp in allExperimenters:
            experimenters.append(exp.tojson())

        allDocumentations = self.documentations()
        documentations = []
        for exp in allDocumentations:
            documentations.append(exp.tojson())

        allDyes = self.dyes()
        dyes = []
        for exp in allDyes:
            dyes.append(exp.tojson())

        allSamples = self.samples()
        samples = []
        for exp in allSamples:
            samples.append(exp.tojson())

        allTargets = self.targets()
        targets = []
        for exp in allTargets:
            targets.append(exp.tojson())

        allTherm_cyc_cons = self.therm_cyc_cons()
        therm_cyc_cons = []
        for exp in allTherm_cyc_cons:
            therm_cyc_cons.append(exp.tojson())

        allExperiments = self.experiments()
        experiments = []
        for exp in allExperiments:
            experiments.append(exp.tojson())

        data = {
            "rdml": {
                "version": self["version"],
                "dateMade": self["dateMade"],
                "dateUpdated": self["dateUpdated"],
                "ids": rdmlids,
                "experimenters": experimenters,
                "documentations": documentations,
                "dyes": dyes,
                "samples": samples,
                "targets": targets,
                "therm_cyc_cons": therm_cyc_cons,
                "experiments": experiments
            }
        }
        return data


class Rdmlid:
    """RDML-Python library

    The rdml id element used to read and edit one experimenter.

    Attributes:
        _node: The rdml id node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an rdml id instance.

        Args:
            self: The class self parameter.
            node: The experimenter node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """

        if key in ["publisher", "serialNumber"]:
            return _get_first_child_text(self._node, key)
        if key in ["MD5Hash"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key in ["publisher", "serialNumber"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["MD5Hash"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["publisher", "serialNumber", "MD5Hash"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "publisher": _get_first_child_text(self._node, "publisher"),
            "serialNumber": _get_first_child_text(self._node, "serialNumber")
        }
        _add_first_child_to_dic(self._node, data, True, "MD5Hash")
        return data


class Experimenter:
    """RDML-Python library

    The experimenter element used to read and edit one experimenter.

    Attributes:
        _node: The experimenter node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an experimenter instance.

        Args:
            self: The class self parameter.
            node: The experimenter node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key in ["firstName", "lastName"]:
            return _get_first_child_text(self._node, key)
        if key in ["email", "labName", "labAddress"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key in ["firstName", "lastName"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["email", "labName", "labAddress"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allTh = _get_all_children(par, "thermalCyclingConditions")
            for node in allTh:
                subNodes = _get_all_children(node, "experimenter")
                for subNode in subNodes:
                    if subNode.attrib['id'] == oldValue:
                        subNode.attrib['id'] = value
            allExp = _get_all_children(par, "experiment")
            for node in allExp:
                subNodes = _get_all_children(node, "run")
                for subNode in subNodes:
                    lastNodes = _get_all_children(subNode, "experimenter")
                    for lastNode in lastNodes:
                        if lastNode.attrib['id'] == oldValue:
                            lastNode.attrib['id'] = value
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "firstName", "lastName", "email", "labName", "labAddress"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
            "firstName": _get_first_child_text(self._node, "firstName"),
            "lastName": _get_first_child_text(self._node, "lastName")
        }
        _add_first_child_to_dic(self._node, data, True, "email")
        _add_first_child_to_dic(self._node, data, True, "labName")
        _add_first_child_to_dic(self._node, data, True, "labAddress")
        return data


class Documentation:
    """RDML-Python library

    The documentation element used to read and edit one documentation tag.

    Attributes:
        _node: The documentation node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an documentation instance.

        Args:
            self: The class self parameter.
            node: The documentation node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the documentation subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "text":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the documentation subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "text":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allSam = _get_all_children(par, "sample")
            for node in allSam:
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    if subNode.attrib['id'] == oldValue:
                        subNode.attrib['id'] = value
            allTh = _get_all_children(par, "target")
            for node in allTh:
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    if subNode.attrib['id'] == oldValue:
                        subNode.attrib['id'] = value
            allTh = _get_all_children(par, "thermalCyclingConditions")
            for node in allTh:
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    if subNode.attrib['id'] == oldValue:
                        subNode.attrib['id'] = value
            allExp = _get_all_children(par, "experiment")
            for node in allExp:
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    if subNode.attrib['id'] == oldValue:
                        subNode.attrib['id'] = value
                subNodes = _get_all_children(node, "run")
                for subNode in subNodes:
                    lastNodes = _get_all_children(subNode, "documentation")
                    for lastNode in lastNodes:
                        if lastNode.attrib['id'] == oldValue:
                            lastNode.attrib['id'] = value
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "text"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "text")
        return data


class Dye:
    """RDML-Python library

    The dye element used to read and edit one dye.

    Attributes:
        _node: The dye node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an dye instance.

        Args:
            self: The class self parameter.
            node: The dye node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the dye subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the dye subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allTar = _get_all_children(par, "target")
            for node in allTar:
                forId = _get_first_child(node, "dyeId")
                if forId is not None:
                    if forId.attrib['id'] == oldValue:
                        forId.attrib['id'] = value
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        return data


class Sample:
    """RDML-Python library

    The samples element used to read and edit one sample.

    Attributes:
        _node: The sample node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an sample instance.

        Args:
            self: The class self parameter.
            node: The sample node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "type":
            return _get_first_child_text(self._node, key)
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key in ["interRunCalibrator", "calibratorSample"]:
            return _get_first_child_bool(self._node, key, triple=True)
        if key in ["cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                   "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions"]:
            ele = _get_first_child(self._node, "cdnaSynthesisMethod")
            if ele is None:
                return None
            if key == "cdnaSynthesisMethod_enzyme":
                return _get_first_child_text(ele, "enzyme")
            if key == "cdnaSynthesisMethod_primingMethod":
                return _get_first_child_text(ele, "primingMethod")
            if key == "cdnaSynthesisMethod_dnaseTreatment":
                return _get_first_child_text(ele, "dnaseTreatment")
            if key == "cdnaSynthesisMethod_thermalCyclingConditions":
                forId = _get_first_child(ele, "thermalCyclingConditions")
                if forId is not None:
                    return forId.attrib['id']
                else:
                    return None
            raise RdmlError('Sample cdnaSynthesisMethod programming read error.')
        if key == "quantity":
            ele = _get_first_child(self._node, key)
            vdic = {}
            vdic["value"] = _get_first_child_text(ele, "value")
            vdic["unit"] = _get_first_child_text(ele, "unit")
            if len(vdic.keys()) != 0:
                return vdic
            else:
                return None
        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            if key in ["templateRNAQuality", "templateDNAQuality"]:
                ele = _get_first_child(self._node, key)
                vdic = {}
                vdic["method"] = _get_first_child_text(ele, "method")
                vdic["result"] = _get_first_child_text(ele, "result")
                if len(vdic.keys()) != 0:
                    return vdic
                else:
                    return None
            if key in ["templateRNAQuantity", "templateDNAQuantity"]:
                ele = _get_first_child(self._node, key)
                vdic = {}
                vdic["value"] = _get_first_child_text(ele, "value")
                vdic["unit"] = _get_first_child_text(ele, "unit")
                if len(vdic.keys()) != 0:
                    return vdic
                else:
                    return None
        if ver == "1.2" or ver == "1.3":
            if key == "templateQuantity":
                ele = _get_first_child(self._node, key)
                vdic = {}
                vdic["nucleotide"] = _get_first_child_text(ele, "nucleotide")
                vdic["conc"] = _get_first_child_text(ele, "conc")
                if len(vdic.keys()) != 0:
                    return vdic
                else:
                    return None
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "type":
            if value not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
                raise RdmlError('Unknown or unsupported sample type value "' + value + '".')

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "type":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key in ["interRunCalibrator", "calibratorSample"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "bool")
        if key in ["cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                   "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions"]:
            ele = _get_or_create_subelement(self._node, "cdnaSynthesisMethod", self.xmlkeys())
            if key == "cdnaSynthesisMethod_enzyme":
                _change_subelement(ele, "enzyme",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "string")
            if key == "cdnaSynthesisMethod_primingMethod":
                if value not in ["", "oligo-dt", "random", "target-specific", "oligo-dt and random", "other"]:
                    raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
                _change_subelement(ele, "primingMethod",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "string")
            if key == "cdnaSynthesisMethod_dnaseTreatment":
                _change_subelement(ele, "dnaseTreatment",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "bool")
            if key == "cdnaSynthesisMethod_thermalCyclingConditions":
                forId = _get_or_create_subelement(ele, "thermalCyclingConditions",
                                                  ["enzyme", "primingMethod", "dnaseTreatment",
                                                   "thermalCyclingConditions"])
                if value is not None and value != "":
                    # We do not check that ID is valid to allow recreate_lost_ids()
                    forId.attrib['id'] = value
                else:
                    ele.remove(forId)
            _remove_irrelevant_subelement(self._node, "cdnaSynthesisMethod")
            return
        if key == "quantity":
            if value is None:
                return
            if "value" not in value or "unit" not in value:
                raise RdmlError('Sample ' + key + ' must have a dictionary with "value" and "unit" as value.')
            if value["unit"] not in ["", "cop", "fold", "dil", "ng", "nMol", "other"]:
                raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
            ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
            _change_subelement(ele, "value", ["value", "unit"], value["value"], True, "float")
            if value["value"] != "":
                _change_subelement(ele, "unit", ["value", "unit"], value["unit"], True, "string")
            else:
                _change_subelement(ele, "unit", ["value", "unit"], "", True, "string")
            _remove_irrelevant_subelement(self._node, key)
            return
        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            if key in ["templateRNAQuality", "templateDNAQuality"]:
                if value is None:
                    return
                if "method" not in value or "result" not in value:
                    raise RdmlError('"' + key + '" must have a dictionary with "method" and "result" as value.')
                ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
                _change_subelement(ele, "method", ["method", "result"], value["method"], True, "string")
                _change_subelement(ele, "result", ["method", "result"], value["result"], True, "float")
                _remove_irrelevant_subelement(self._node, key)
                return
            if key in ["templateRNAQuantity", "templateDNAQuantity"]:
                if value is None:
                    return
                if "value" not in value or "unit" not in value:
                    raise RdmlError('Sample ' + key + ' must have a dictionary with "value" and "unit" as value.')
                if value["unit"] not in ["", "cop", "fold", "dil", "ng", "nMol", "other"]:
                    raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
                ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
                _change_subelement(ele, "value", ["value", "unit"], value["value"], True, "float")
                if value["value"] != "":
                    _change_subelement(ele, "unit", ["value", "unit"], value["unit"], True, "string")
                else:
                    _change_subelement(ele, "unit", ["value", "unit"], "", True, "string")
                _remove_irrelevant_subelement(self._node, key)
                return
        if ver == "1.2" or ver == "1.3":
            if key == "templateQuantity":
                if value is None:
                    return
                if "nucleotide" not in value or "conc" not in value:
                    raise RdmlError('Sample ' + key + ' must have a dictionary with "nucleotide" and "conc" as value.')
                if value["nucleotide"] not in ["", "DNA", "genomic DNA", "cDNA", "RNA"]:
                    raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
                ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
                _change_subelement(ele, "conc", ["conc", "nucleotide"], value["conc"], True, "float")
                if value["conc"] != "":
                    _change_subelement(ele, "nucleotide", ["conc", "nucleotide"], value["nucleotide"], True, "string")
                else:
                    _change_subelement(ele, "nucleotide", ["conc", "nucleotide"], "", True, "string")
                _remove_irrelevant_subelement(self._node, key)
                return
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allExp = _get_all_children(par, "experiment")
            for node in allExp:
                subNodes = _get_all_children(node, "run")
                for subNode in subNodes:
                    reactNodes = _get_all_children(subNode, "react")
                    for reactNode in reactNodes:
                        lastNodes = _get_all_children(reactNode, "sample")
                        for lastNode in lastNodes:
                            if lastNode.attrib['id'] == oldValue:
                                lastNode.attrib['id'] = value
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return ["id", "description", "type", "interRunCalibrator", "quantity", "calibratorSample",
                    "cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                    "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions",
                    "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]
        return ["id", "description", "annotation", "type", "interRunCalibrator", "quantity", "calibratorSample",
                "cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions",
                "templateQuantity"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return ["description", "documentation", "xRef", "type", "interRunCalibrator",
                    "quantity", "calibratorSample", "cdnaSynthesisMethod",
                    "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]
        return ["description", "documentation", "xRef", "annotation", "type", "interRunCalibrator",
                "quantity", "calibratorSample", "cdnaSynthesisMethod", "templateQuantity"]

    def xrefs(self):
        """Returns a list of the xrefs in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with name and id strings.
        """

        xref = _get_all_children(self._node, "xRef")
        ret = []
        for node in xref:
            data = {}
            _add_first_child_to_dic(node, data, True, "name")
            _add_first_child_to_dic(node, data, True, "id")
            ret.append(data)
        return ret

    def new_xref(self, name=None, id=None, newposition=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            name: Publisher who created the xRef
            id: Serial Number for this sample provided by publisher
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if name is None and id is None:
            raise RdmlError('Either name or id is required to create a xRef.')
        new_node = ET.Element("xRef")
        _add_new_subelement(new_node, "xRef", "name", name, True)
        _add_new_subelement(new_node, "xRef", "id", id, True)
        place = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_xref(self, oldposition, newposition=None, name=None, id=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            name: Publisher who created the xRef
            id: Serial Number for this sample provided by publisher

        Returns:
            Nothing, changes self.
        """

        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a xRef.')
        if (name is None or name == "") and (id is None or id == ""):
            self.delete_xref(oldposition)
            return
        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        _change_subelement(ele, "name", ["name", "id"], name, True, "string")
        _change_subelement(ele, "id", ["name", "id"], id, True, "string", id_as_element=True)
        self._node.insert(pos, ele)

    def move_xref(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        self._node.insert(pos, ele)

    def delete_xref(self, byposition):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "xRef", None, byposition)
        self._node.remove(elem)

    def annotations(self):
        """Returns a list of the annotations in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with property and value strings.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return []
        xref = _get_all_children(self._node, "annotation")
        ret = []
        for node in xref:
            data = {}
            _add_first_child_to_dic(node, data, True, "property")
            _add_first_child_to_dic(node, data, True, "value")
            ret.append(data)
        return ret

    def new_annotation(self, property=None, value=None, newposition=None):
        """Creates a new annotation element.

        Args:
            self: The class self parameter.
            property: The property
            value: Its value
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return
        if property is None or value is None:
            raise RdmlError('Property and value are required to create a annotation.')
        new_node = ET.Element("annotation")
        _add_new_subelement(new_node, "annotation", "property", property, True)
        _add_new_subelement(new_node, "annotation", "value", value, True)
        place = _get_tag_pos(self._node, "annotation", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_annotation(self, oldposition, newposition=None, property=None, value=None):
        """Edits an annotation element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            property: The property
            value: Its value

        Returns:
            Nothing, changes self.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return
        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a annotation.')
        if (property is None or property == "") or (value is None or value == ""):
            self.delete_annotation(oldposition)
            return
        pos = _get_tag_pos(self._node, "annotation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "annotation", None, oldposition)
        _change_subelement(ele, "property", ["property", "value"], property, True, "string")
        _change_subelement(ele, "value", ["property", "value"], value, True, "string")
        self._node.insert(pos, ele)

    def move_annotation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return
        pos = _get_tag_pos(self._node, "annotation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "annotation", None, oldposition)
        self._node.insert(pos, ele)

    def delete_annotation(self, byposition):
        """Deletes an annotation element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return
        elem = _get_first_child_by_pos_or_id(self._node, "annotation", None, byposition)
        self._node.remove(elem)

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """
        par = self._node.getparent()
        ver = par.get('version')

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["xRefs"] = self.xrefs()
        if ver == "1.2" or ver == "1.3":
            data["annotations"] = self.annotations()
        _add_first_child_to_dic(self._node, data, False, "type")
        _add_first_child_to_dic(self._node, data, True, "interRunCalibrator")
        elem = _get_first_child(self._node, "quantity")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "value")
            _add_first_child_to_dic(elem, qdic, False, "unit")
            data["quantity"] = qdic
        _add_first_child_to_dic(self._node, data, True, "calibratorSample")
        elem = _get_first_child(self._node, "cdnaSynthesisMethod")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, True, "enzyme")
            _add_first_child_to_dic(elem, qdic, True, "primingMethod")
            _add_first_child_to_dic(elem, qdic, True, "dnaseTreatment")
            forId = _get_first_child(elem, "thermalCyclingConditions")
            if forId is not None:
                if forId.attrib['id'] != "":
                    qdic["thermalCyclingConditions"] = forId.attrib['id']
            if len(qdic.keys()) != 0:
                data["cdnaSynthesisMethod"] = qdic
        if ver == "1.1":
            elem = _get_first_child(self._node, "templateRNAQuantity")
            if elem is not None:
                qdic = {}
                _add_first_child_to_dic(elem, qdic, False, "value")
                _add_first_child_to_dic(elem, qdic, False, "unit")
                data["templateRNAQuantity"] = qdic
            elem = _get_first_child(self._node, "templateRNAQuality")
            if elem is not None:
                qdic = {}
                _add_first_child_to_dic(elem, qdic, False, "method")
                _add_first_child_to_dic(elem, qdic, False, "result")
                data["templateRNAQuality"] = qdic
            elem = _get_first_child(self._node, "templateDNAQuantity")
            if elem is not None:
                qdic = {}
                _add_first_child_to_dic(elem, qdic, False, "value")
                _add_first_child_to_dic(elem, qdic, False, "unit")
                data["templateDNAQuantity"] = qdic
            elem = _get_first_child(self._node, "templateDNAQuality")
            if elem is not None:
                qdic = {}
                _add_first_child_to_dic(elem, qdic, False, "method")
                _add_first_child_to_dic(elem, qdic, False, "result")
                data["templateDNAQuality"] = qdic
        if ver == "1.2" or ver == "1.3":
            elem = _get_first_child(self._node, "templateQuantity")
            if elem is not None:
                qdic = {}
                _add_first_child_to_dic(elem, qdic, False, "nucleotide")
                _add_first_child_to_dic(elem, qdic, False, "conc")
                data["templateQuantity"] = qdic
        return data


class Target:
    """RDML-Python library

    The target element used to read and edit one target.

    Attributes:
        _node: The target node of the RDML XML object.
        _rdmlFilename: The RDML filename
    """

    def __init__(self, node, rdmlFilename):
        """Inits an target instance.

        Args:
            self: The class self parameter.
            node: The target node.
            rdmlFilename: The RDML filename.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlFilename = rdmlFilename

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the target subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "type":
            return _get_first_child_text(self._node, key)
        if key in ["description", "amplificationEfficiencyMethod", "amplificationEfficiency", "detectionLimit"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key == "dyeId":
            forId = _get_first_child(self._node, key)
            if forId is not None:
                return forId.attrib['id']
            else:
                return None
        if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                   "sequences_forwardPrimer_sequence", "sequences_reversePrimer_threePrimeTag",
                   "sequences_reversePrimer_fivePrimeTag", "sequences_reversePrimer_sequence",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                   "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence"]:
            prim = _get_first_child(self._node, "sequences")
            if prim is None:
                return None
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                sec = _get_first_child(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                sec = _get_first_child(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                sec = _get_first_child(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                sec = _get_first_child(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                sec = _get_first_child(prim, "amplicon")
            if sec is None:
                return None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_reversePrimer_threePrimeTag",
                       "sequences_probe1_threePrimeTag", "sequences_probe2_threePrimeTag",
                       "sequences_amplicon_threePrimeTag"]:
                return _get_first_child_text(sec, "threePrimeTag")
            if key in ["sequences_forwardPrimer_fivePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_probe1_fivePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_amplicon_fivePrimeTag"]:
                return _get_first_child_text(sec, "fivePrimeTag")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                return _get_first_child_text(sec, "sequence")
            raise RdmlError('Target sequences programming read error.')
        if key in ["commercialAssay_company", "commercialAssay_orderNumber"]:
            prim = _get_first_child(self._node, "commercialAssay")
            if prim is None:
                return None
            if key == "commercialAssay_company":
                return _get_first_child_text(prim, "company")
            if key == "commercialAssay_orderNumber":
                return _get_first_child_text(prim, "orderNumber")
        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.2" or ver == "1.3":
            if key == "amplificationEfficiencySE":
                var = _get_first_child_text(self._node, key)
                if var == "":
                    return None
                else:
                    return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the target subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "type":
            if value not in ["ref", "toi"]:
                raise RdmlError('Unknown or unsupported target type value "' + value + '".')

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "type":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["description", "amplificationEfficiencyMethod"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key in ["amplificationEfficiency", "detectionLimit"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        if key == "dyeId":
            forId = _get_or_create_subelement(self._node, "dyeId", self.xmlkeys())
            if value is not None and value != "":
                # We do not check that ID is valid to allow recreate_lost_ids()
                forId.attrib['id'] = value
            else:
                self._node.remove(forId)
            return
        if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                   "sequences_forwardPrimer_sequence", "sequences_reversePrimer_threePrimeTag",
                   "sequences_reversePrimer_fivePrimeTag", "sequences_reversePrimer_sequence",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                   "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence"]:
            prim = _get_or_create_subelement(self._node, "sequences", self.xmlkeys())
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                sec = _get_or_create_subelement(prim, "forwardPrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                sec = _get_or_create_subelement(prim, "reversePrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                sec = _get_or_create_subelement(prim, "probe1",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                sec = _get_or_create_subelement(prim, "probe2",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                sec = _get_or_create_subelement(prim, "amplicon",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if sec is None:
                return None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_reversePrimer_threePrimeTag",
                       "sequences_probe1_threePrimeTag", "sequences_probe2_threePrimeTag",
                       "sequences_amplicon_threePrimeTag"]:
                _change_subelement(sec, "threePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_fivePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_probe1_fivePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_amplicon_fivePrimeTag"]:
                _change_subelement(sec, "fivePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                _change_subelement(sec, "sequence",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                _remove_irrelevant_subelement(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                _remove_irrelevant_subelement(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                _remove_irrelevant_subelement(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                _remove_irrelevant_subelement(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                _remove_irrelevant_subelement(prim, "amplicon")
            _remove_irrelevant_subelement(self._node, "sequences")
            return
        if key in ["commercialAssay_company", "commercialAssay_orderNumber"]:
            ele = _get_or_create_subelement(self._node, "commercialAssay", self.xmlkeys())
            if key == "commercialAssay_company":
                _change_subelement(ele, "company", ["company", "orderNumber"], value, True, "string")
            if key == "commercialAssay_orderNumber":
                _change_subelement(ele, "orderNumber", ["company", "orderNumber"], value, True, "string")
            _remove_irrelevant_subelement(self._node, "commercialAssay")
            return
        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.2" or ver == "1.3":
            if key == "amplificationEfficiencySE":
                return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allExp = _get_all_children(par, "experiment")
            for node in allExp:
                subNodes = _get_all_children(node, "run")
                for subNode in subNodes:
                    reactNodes = _get_all_children(subNode, "react")
                    for reactNode in reactNodes:
                        dataNodes = _get_all_children(reactNode, "data")
                        for dataNode in dataNodes:
                            lastNodes = _get_all_children(dataNode, "tar")
                            for lastNode in lastNodes:
                                if lastNode.attrib['id'] == oldValue:
                                    lastNode.attrib['id'] = value
                        partit = _get_first_child(reactNode, "partitions")
                        if partit is not None:
                            digDataNodes = _get_all_children(partit, "data")
                            for digDataNode in digDataNodes:
                                lastNodes = _get_all_children(digDataNode, "tar")
                                for lastNode in lastNodes:
                                    if lastNode.attrib['id'] == oldValue:
                                        lastNode.attrib['id'] = value

            # Search in Table files
            if self._rdmlFilename is not None and self._rdmlFilename != "":
                if zipfile.is_zipfile(self._rdmlFilename):
                    fileList = []
                    tempName = ""
                    flipFiles = False
                    with zipfile.ZipFile(self._rdmlFilename, 'r') as RDMLin:
                        for item in RDMLin.infolist():
                            if re.search("^partitions/", item.filename):
                                fileContent = RDMLin.read(item.filename).decode('utf-8')
                                newlineFix = fileContent.replace("\r\n", "\n")
                                tabLines = newlineFix.split("\n")
                                header = tabLines[0].split("\t")
                                needRewrite = False
                                for cell in header:
                                    if cell == oldValue:
                                        needRewrite = True
                                if needRewrite:
                                    fileList.append(item.filename)
                        if len(fileList) > 0:
                            tempFolder, tempName = tempfile.mkstemp(dir=os.path.dirname(self._rdmlFilename))
                            os.close(tempFolder)
                            flipFiles = True
                            with zipfile.ZipFile(tempName, mode='w', compression=zipfile.ZIP_DEFLATED) as RDMLout:
                                RDMLout.comment = RDMLin.comment
                                for item in RDMLin.infolist():
                                    if item.filename not in fileList:
                                        RDMLout.writestr(item, RDMLin.read(item.filename))
                                    else:
                                        fileContent = RDMLin.read(item.filename).decode('utf-8')
                                        newlineFix = fileContent.replace("\r\n", "\n")
                                        tabLines = newlineFix.split("\n")
                                        header = tabLines[0].split("\t")
                                        headerText = ""
                                        for cell in header:
                                            if cell == oldValue:
                                                headerText += value + "\t"
                                            else:
                                                headerText += cell + "\t"
                                        outFileStr = re.sub(r'\t$', '\n', headerText)
                                        for tabLine in tabLines[1:]:
                                            if tabLine != "":
                                                outFileStr += tabLine + "\n"
                                        RDMLout.writestr(item.filename, outFileStr)
                    if flipFiles:
                        os.remove(self._rdmlFilename)
                        os.rename(tempName, self._rdmlFilename)
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "type", "amplificationEfficiencyMethod", "amplificationEfficiency",
                "detectionLimit", "dyeId", "sequences_forwardPrimer_threePrimeTag",
                "sequences_forwardPrimer_fivePrimeTag", "sequences_forwardPrimer_sequence",
                "sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                "sequences_reversePrimer_sequence", "sequences_probe1_threePrimeTag",
                "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence", "sequences_amplicon_threePrimeTag",
                "sequences_amplicon_fivePrimeTag", "sequences_amplicon_sequence", "commercialAssay_company",
                "commercialAssay_orderNumber"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "xRef", "type", "amplificationEfficiencyMethod",
                "amplificationEfficiency", "detectionLimit", "dyeId", "sequences", "commercialAssay"]

    def xrefs(self):
        """Returns a list of the xrefs in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with name and id strings.
        """

        xref = _get_all_children(self._node, "xRef")
        ret = []
        for node in xref:
            data = {}
            _add_first_child_to_dic(node, data, True, "name")
            _add_first_child_to_dic(node, data, True, "id")
            ret.append(data)
        return ret

    def new_xref(self, name=None, id=None, newposition=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            name: Publisher who created the xRef
            id: Serial Number for this target provided by publisher
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if name is None and id is None:
            raise RdmlError('Either name or id is required to create a xRef.')
        new_node = ET.Element("xRef")
        _add_new_subelement(new_node, "xRef", "name", name, True)
        _add_new_subelement(new_node, "xRef", "id", id, True)
        place = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_xref(self, oldposition, newposition=None, name=None, id=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            name: Publisher who created the xRef
            id: Serial Number for this target provided by publisher

        Returns:
            Nothing, changes self.
        """

        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a xRef.')
        if (name is None or name == "") and (id is None or id == ""):
            self.delete_xref(oldposition)
            return
        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        _change_subelement(ele, "name", ["name", "id"], name, True, "string")
        _change_subelement(ele, "id", ["name", "id"], id, True, "string", id_as_element=True)
        self._node.insert(pos, ele)

    def move_xref(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        self._node.insert(pos, ele)

    def delete_xref(self, byposition):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "xRef", None, byposition)
        self._node.remove(elem)

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["xRefs"] = self.xrefs()
        _add_first_child_to_dic(self._node, data, False, "type")
        _add_first_child_to_dic(self._node, data, True, "amplificationEfficiencyMethod")
        _add_first_child_to_dic(self._node, data, True, "amplificationEfficiency")
        _add_first_child_to_dic(self._node, data, True, "detectionLimit")
        forId = _get_first_child(self._node, "dyeId")
        if forId is not None:
            if forId.attrib['id'] != "":
                data["dyeId"] = forId.attrib['id']
        elem = _get_first_child(self._node, "sequences")
        if elem is not None:
            qdic = {}
            sec = _get_first_child(elem, "forwardPrimer")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["forwardPrimer"] = sdic
            sec = _get_first_child(elem, "reversePrimer")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["reversePrimer"] = sdic
            sec = _get_first_child(elem, "probe1")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["probe1"] = sdic
            sec = _get_first_child(elem, "probe2")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["probe2"] = sdic
            sec = _get_first_child(elem, "amplicon")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["amplicon"] = sdic
            if len(qdic.keys()) != 0:
                data["sequences"] = qdic
        elem = _get_first_child(self._node, "commercialAssay")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, True, "company")
            _add_first_child_to_dic(elem, qdic, True, "orderNumber")
            if len(qdic.keys()) != 0:
                data["commercialAssay"] = qdic
        return data


class Therm_cyc_cons:
    """RDML-Python library

    The thermalCyclingConditions element used to read and edit one thermal Cycling Conditions.

    Attributes:
        _node: The thermalCyclingConditions node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an thermalCyclingConditions instance.

        Args:
            self: The class self parameter.
            node: The thermalCyclingConditions node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the thermalCyclingConditions subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key in ["description", "lidTemperature"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the thermalCyclingConditions subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key == "lidTemperature":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        raise KeyError

    def change_id(self, value, merge_with_id=False):
        """Changes the value for the id.

        Args:
            self: The class self parameter.
            value: The new value for the id.
            merge_with_id: If True only allow a unique id, if False only rename its uses with existing id.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        oldValue = self._node.get('id')
        if oldValue != value:
            par = self._node.getparent()
            if not _string_to_bool(merge_with_id, triple=False):
                _change_subelement(self._node, "id", self.xmlkeys(), value, False, "string")
            else:
                groupTag = self._node.tag.replace("{http://www.rdml.org}", "")
                if _check_unique_id(par, groupTag, value):
                    raise RdmlError('The ' + groupTag + ' id "' + value + '" does not exist.')
            allSam = _get_all_children(par, "sample")
            for node in allSam:
                subNode = _get_first_child(node, "cdnaSynthesisMethod")
                if subNode is not None:
                    forId = _get_first_child(subNode, "thermalCyclingConditions")
                    if forId is not None:
                        if forId.attrib['id'] == oldValue:
                            forId.attrib['id'] = value
            allExp = _get_all_children(par, "experiment")
            for node in allExp:
                subNodes = _get_all_children(node, "run")
                for subNode in subNodes:
                    forId = _get_first_child(subNode, "thermalCyclingConditions")
                    if forId is not None:
                        if forId.attrib['id'] == oldValue:
                            forId.attrib['id'] = value
        return

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "lidTemperature"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "lidTemperature", "experimenter", "step"]

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def experimenter_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "experimenter")

    def update_experimenter_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.experimenter_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "experimenter", id)
                    place = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "experimenter", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_experimenter(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "experimenter", None, oldposition)
        self._node.insert(pos, ele)

    def steps(self):
        """Returns a list of all step elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all step elements.
        """

        # The steps are sorted transiently to not modify the file in a read situation
        exp = _get_all_children(self._node, "step")
        srt_exp = sorted(exp, key=_get_step_sort_nr)
        ret = []
        for node in srt_exp:
            ret.append(Step(node))
        return ret

    def new_step_temperature(self, temperature, duration,
                             temperatureChange=None, durationChange=None,
                             measure=None, ramp=None, nr=None):
        """Creates a new step element.

        Args:
            self: The class self parameter.
            temperature: The temperature of the step in degrees Celsius (required)
            duration: The duration of this step in seconds (required)
            temperatureChange: The change of the temperature from one cycle to the next (optional)
            durationChange: The change of the duration from one cycle to the next (optional)
            measure: Indicates to make a measurement and store it as meltcurve or real-time data (optional)
            ramp: Limit temperature change from one step to the next in degrees Celsius per second (optional)
            nr: Step unique nr (optional)

        Returns:
            Nothing, changes self.
        """

        if measure is not None and measure not in ["", "real time", "meltcurve"]:
            raise RdmlError('Unknown or unsupported step measure value: "' + measure + '".')
        nr = int(nr)
        count = _get_number_of_children(self._node, "step")
        new_node = ET.Element("step")
        xml_temp_step = ["temperature", "duration", "temperatureChange", "durationChange", "measure", "ramp"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = ET.SubElement(new_node, "temperature")
        _change_subelement(subel, "temperature", xml_temp_step, temperature, False, "float")
        _change_subelement(subel, "duration", xml_temp_step, duration, False, "posint")
        _change_subelement(subel, "temperatureChange", xml_temp_step, temperatureChange, True, "float")
        _change_subelement(subel, "durationChange", xml_temp_step, durationChange, True, "int")
        _change_subelement(subel, "measure", xml_temp_step, measure, True, "string")
        _change_subelement(subel, "ramp", xml_temp_step, ramp, True, "float")
        place = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + count
        self._node.insert(place, new_node)
        # Now move step at final position
        self.move_step(count + 1, nr)

    def new_step_gradient(self, highTemperature, lowTemperature, duration,
                          temperatureChange=None, durationChange=None,
                          measure=None, ramp=None, nr=None):
        """Creates a new step element.

        Args:
            self: The class self parameter.
            highTemperature: The high gradient temperature of the step in degrees Celsius (required)
            lowTemperature: The low gradient temperature of the step in degrees Celsius (required)
            duration: The duration of this step in seconds (required)
            temperatureChange: The change of the temperature from one cycle to the next (optional)
            durationChange: The change of the duration from one cycle to the next (optional)
            measure: Indicates to make a measurement and store it as meltcurve or real-time data (optional)
            ramp: Limit temperature change from one step to the next in degrees Celsius per second (optional)
            nr: Step unique nr (optional)

        Returns:
            Nothing, changes self.
        """

        if measure is not None and measure not in ["", "real time", "meltcurve"]:
            raise RdmlError('Unknown or unsupported step measure value: "' + measure + '".')
        nr = int(nr)
        count = _get_number_of_children(self._node, "step")
        new_node = ET.Element("step")
        xml_temp_step = ["highTemperature", "lowTemperature", "duration", "temperatureChange",
                         "durationChange", "measure", "ramp"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = ET.SubElement(new_node, "gradient")
        _change_subelement(subel, "highTemperature", xml_temp_step, highTemperature, False, "float")
        _change_subelement(subel, "lowTemperature", xml_temp_step, lowTemperature, False, "float")
        _change_subelement(subel, "duration", xml_temp_step, duration, False, "posint")
        _change_subelement(subel, "temperatureChange", xml_temp_step, temperatureChange, True, "float")
        _change_subelement(subel, "durationChange", xml_temp_step, durationChange, True, "int")
        _change_subelement(subel, "measure", xml_temp_step, measure, True, "string")
        _change_subelement(subel, "ramp", xml_temp_step, ramp, True, "float")
        place = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + count
        self._node.insert(place, new_node)
        # Now move step at final position
        self.move_step(count + 1, nr)

    def new_step_loop(self, goto, repeat, nr=None):
        """Creates a new step element.

        Args:
            self: The class self parameter.
            goto: The step nr to go back to (required)
            repeat: The number of times to go back to goto step, one less than cycles (optional)
            nr: Step unique nr (optional)

        Returns:
            Nothing, changes self.
        """

        nr = int(nr)
        count = _get_number_of_children(self._node, "step")
        new_node = ET.Element("step")
        xml_temp_step = ["goto", "repeat"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = ET.SubElement(new_node, "loop")
        _change_subelement(subel, "goto", xml_temp_step, goto, False, "posint")
        _change_subelement(subel, "repeat", xml_temp_step, repeat, False, "posint")
        place = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + count
        self._node.insert(place, new_node)
        # Now move step at final position
        self.move_step(count + 1, nr)

    def new_step_pause(self, temperature, nr=None):
        """Creates a new step element.

        Args:
            self: The class self parameter.
            temperature: The temperature of the step in degrees Celsius (required)
            nr: Step unique nr (optional)

        Returns:
            Nothing, changes self.
        """

        nr = int(nr)
        count = _get_number_of_children(self._node, "step")
        new_node = ET.Element("step")
        xml_temp_step = ["temperature"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = ET.SubElement(new_node, "pause")
        _change_subelement(subel, "temperature", xml_temp_step, temperature, False, "float")
        place = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + count
        self._node.insert(place, new_node)
        # Now move step at final position
        self.move_step(count + 1, nr)

    def new_step_lidOpen(self, nr=None):
        """Creates a new step element.

        Args:
            self: The class self parameter.
            nr: Step unique nr (optional)

        Returns:
            Nothing, changes self.
        """

        nr = int(nr)
        count = _get_number_of_children(self._node, "step")
        new_node = ET.Element("step")
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        ET.SubElement(new_node, "lidOpen")
        place = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + count
        self._node.insert(place, new_node)
        # Now move step at final position
        self.move_step(count + 1, nr)

    def cleanup_steps(self):
        """The steps may not be in a order that makes sense. This function fixes it.

        Args:
            self: The class self parameter.

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        # The steps in the xml may be not sorted by "nr", so sort first
        exp = _get_all_children(self._node, "step")
        srt_exp = sorted(exp, key=_get_step_sort_nr)
        i = 0
        for node in srt_exp:
            if _get_step_sort_nr(node) != _get_step_sort_nr(exp[i]):
                pos = _get_first_tag_pos(self._node, "step", self.xmlkeys()) + i
                self._node.insert(pos, node)
            i += 1

        # The steps in the xml may not have the correct numbering, so fix it
        exp = _get_all_children(self._node, "step")
        i = 1
        for node in exp:
            if _get_step_sort_nr(node) != i:
                elem = _get_first_child(node, "nr")
                elem.text = str(i)
            i += 1

    def move_step(self, oldnr, newnr):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldnr: The old position of the element
            newnr: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        # The steps in the xml may be not sorted well, so fix it
        self.cleanup_steps()

        # Change the nr
        _move_subelement_pos(self._node, "step", oldnr - 1, self.xmlkeys(), newnr - 1)

        # Fix the nr
        exp = _get_all_children(self._node, "step")
        i = 1
        goto_mod = 0
        goto_start = newnr
        goto_end = oldnr
        if oldnr > newnr:
            goto_mod = 1
        if oldnr < newnr:
            goto_mod = -1
            goto_start = oldnr
            goto_end = newnr
        for node in exp:
            if _get_step_sort_nr(node) != i:
                elem = _get_first_child(node, "nr")
                elem.text = str(i)
            # Fix the goto steps
            ele_type = _get_first_child(node, "loop")
            if ele_type is not None:
                ele_goto = _get_first_child(ele_type, "goto")
                if ele_goto is not None:
                    jump_to = int(ele_goto.text)
                    if goto_start <= jump_to < goto_end:
                        ele_goto.text = str(jump_to + goto_mod)
            i += 1

    def get_step(self, bystep):
        """Returns an sample element by position or id.

        Args:
            self: The class self parameter.
            bystep: Select the element by step nr in the list.

        Returns:
            The found element or None.
        """

        return Step(_get_first_child_by_pos_or_id(self._node, "step", None, bystep - 1))

    def delete_step(self, bystep=None):
        """Deletes an step element.

        Args:
            self: The class self parameter.
            bystep: Select the element by step nr in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "step", None, bystep - 1)
        self._node.remove(elem)
        self.cleanup_steps()
        # Fix the goto steps
        exp = _get_all_children(self._node, "step")
        for node in exp:
            ele_type = _get_first_child(node, "loop")
            if ele_type is not None:
                ele_goto = _get_first_child(ele_type, "goto")
                if ele_goto is not None:
                    jump_to = int(ele_goto.text)
                    if bystep < jump_to:
                        ele_goto.text = str(jump_to - 1)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        allSteps = self.steps()
        steps = []
        for exp in allSteps:
            steps.append(exp.tojson())

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        _add_first_child_to_dic(self._node, data, True, "lidTemperature")
        data["experimenters"] = self.experimenter_ids()
        data["steps"] = steps
        return data


class Step:
    """RDML-Python library

    The samples element used to read and edit one sample.

    Attributes:
        _node: The sample node of the RDML XML object.
    """

    def __init__(self, node):
        """Inits an sample instance.

        Args:
            self: The class self parameter.
            node: The sample node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement. Be aware that change of type deletes all entries
                 except nr and description

        Returns:
            A string of the data or None.
        """

        if key == "nr":
            return _get_first_child_text(self._node, key)
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        ele_type = _get_first_child(self._node, "temperature")
        if ele_type is not None:
            if key == "type":
                return "temperature"
            if key in ["temperature", "duration"]:
                return _get_first_child_text(ele_type, key)
            if key in ["temperatureChange", "durationChange", "measure", "ramp"]:
                var = _get_first_child_text(ele_type, key)
                if var == "":
                    return None
                else:
                    return var
        ele_type = _get_first_child(self._node, "gradient")
        if ele_type is not None:
            if key == "type":
                return "gradient"
            if key in ["highTemperature", "lowTemperature", "duration"]:
                return _get_first_child_text(ele_type, key)
            if key in ["temperatureChange", "durationChange", "measure", "ramp"]:
                var = _get_first_child_text(ele_type, key)
                if var == "":
                    return None
                else:
                    return var
        ele_type = _get_first_child(self._node, "loop")
        if ele_type is not None:
            if key == "type":
                return "loop"
            if key in ["goto", "repeat"]:
                return _get_first_child_text(ele_type, key)
        ele_type = _get_first_child(self._node, "pause")
        if ele_type is not None:
            if key == "type":
                return "pause"
            if key == "temperature":
                return _get_first_child_text(ele_type, key)
        ele_type = _get_first_child(self._node, "lidOpen")
        if ele_type is not None:
            if key == "type":
                return "lidOpen"
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key in ["nr", "type"]:
            raise RdmlError('"' + key + '" can not be set. Use thermal cycling conditions methods instead')
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        ele_type = _get_first_child(self._node, "temperature")
        if ele_type is not None:
            xml_temp_step = ["temperature", "duration", "temperatureChange", "durationChange", "measure", "ramp"]
            if key == "temperature":
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "float")
            if key == "duration":
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "posint")
            if key in ["temperatureChange", "ramp"]:
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "float")
            if key == "durationChange":
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "int")
            if key == "measure":
                if value not in ["", "real time", "meltcurve"]:
                    raise RdmlError('Unknown or unsupported step measure value: "' + value + '".')
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "string")
        ele_type = _get_first_child(self._node, "gradient")
        if ele_type is not None:
            xml_temp_step = ["highTemperature", "lowTemperature", "duration", "temperatureChange",
                             "durationChange", "measure", "ramp"]
            if key in ["highTemperature", "lowTemperature"]:
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "float")
            if key == "duration":
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "posint")
            if key in ["temperatureChange", "ramp"]:
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "float")
            if key == "durationChange":
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "int")
            if key == "measure":
                if value not in ["real time", "meltcurve"]:
                    raise RdmlError('Unknown or unsupported step measure value: "' + value + '".')
                return _change_subelement(ele_type, key, xml_temp_step, value, True, "string")
        ele_type = _get_first_child(self._node, "loop")
        if ele_type is not None:
            xml_temp_step = ["goto", "repeat"]
            if key in xml_temp_step:
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "posint")
        ele_type = _get_first_child(self._node, "pause")
        if ele_type is not None:
            xml_temp_step = ["temperature"]
            if key == "temperature":
                return _change_subelement(ele_type, key, xml_temp_step, value, False, "float")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        ele_type = _get_first_child(self._node, "temperature")
        if ele_type is not None:
            return ["nr", "type", "description", "temperature", "duration", "temperatureChange",
                    "durationChange", "measure", "ramp"]
        ele_type = _get_first_child(self._node, "gradient")
        if ele_type is not None:
            return ["nr", "type", "description", "highTemperature", "lowTemperature", "duration",
                    "temperatureChange", "durationChange", "measure", "ramp"]
        ele_type = _get_first_child(self._node, "loop")
        if ele_type is not None:
            return ["nr", "type", "description", "goto", "repeat"]
        ele_type = _get_first_child(self._node, "pause")
        if ele_type is not None:
            return ["nr", "type", "description", "temperature"]
        ele_type = _get_first_child(self._node, "lidOpen")
        if ele_type is not None:
            return ["nr", "type", "description"]
        return []

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        ele_type = _get_first_child(self._node, "temperature")
        if ele_type is not None:
            return ["temperature", "duration", "temperatureChange", "durationChange", "measure", "ramp"]
        ele_type = _get_first_child(self._node, "gradient")
        if ele_type is not None:
            return ["highTemperature", "lowTemperature", "duration", "temperatureChange",
                    "durationChange", "measure", "ramp"]
        ele_type = _get_first_child(self._node, "loop")
        if ele_type is not None:
            return ["goto", "repeat"]
        ele_type = _get_first_child(self._node, "pause")
        if ele_type is not None:
            return ["temperature"]
        ele_type = _get_first_child(self._node, "lidOpen")
        if ele_type is not None:
            return []
        return []

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {}
        _add_first_child_to_dic(self._node, data, False, "nr")
        _add_first_child_to_dic(self._node, data, True, "description")
        elem = _get_first_child(self._node, "temperature")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "temperature")
            _add_first_child_to_dic(elem, qdic, False, "duration")
            _add_first_child_to_dic(elem, qdic, True, "temperatureChange")
            _add_first_child_to_dic(elem, qdic, True, "durationChange")
            _add_first_child_to_dic(elem, qdic, True, "measure")
            _add_first_child_to_dic(elem, qdic, True, "ramp")
            data["temperature"] = qdic
        elem = _get_first_child(self._node, "gradient")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "highTemperature")
            _add_first_child_to_dic(elem, qdic, False, "lowTemperature")
            _add_first_child_to_dic(elem, qdic, False, "duration")
            _add_first_child_to_dic(elem, qdic, True, "temperatureChange")
            _add_first_child_to_dic(elem, qdic, True, "durationChange")
            _add_first_child_to_dic(elem, qdic, True, "measure")
            _add_first_child_to_dic(elem, qdic, True, "ramp")
            data["gradient"] = qdic
        elem = _get_first_child(self._node, "loop")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "goto")
            _add_first_child_to_dic(elem, qdic, False, "repeat")
            data["loop"] = qdic
        elem = _get_first_child(self._node, "pause")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "temperature")
            data["pause"] = qdic
        elem = _get_first_child(self._node, "lidOpen")
        if elem is not None:
            data["lidOpen"] = "lidOpen"
        return data


class Experiment:
    """RDML-Python library

    The target element used to read and edit one experiment.

    Attributes:
        _node: The target node of the RDML XML object.
        _rdmlFilename: The RDML filename
    """

    def __init__(self, node, rdmlFilename):
        """Inits an experiment instance.

        Args:
            self: The class self parameter.
            node: The experiment node.
            rdmlFilename: The RDML filename.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlFilename = rdmlFilename

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experiment subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the target subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "id":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "run"]

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def runs(self):
        """Returns a list of all run elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all run elements.
        """

        exp = _get_all_children(self._node, "run")
        ret = []
        for node in exp:
            ret.append(Run(node, self._rdmlFilename))
        return ret

    def new_run(self, id, newposition=None):
        """Creates a new run element.

        Args:
            self: The class self parameter.
            id: Run unique id (required)
            newposition: Run position in the list of experiments (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = _create_new_element(self._node, "run", id)
        place = _get_tag_pos(self._node, "run", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_run(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Run unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "run", id, self.xmlkeys(), newposition)

    def get_run(self, byid=None, byposition=None):
        """Returns an run element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Run(_get_first_child_by_pos_or_id(self._node, "run", byid, byposition), self._rdmlFilename)

    def delete_run(self, byid=None, byposition=None):
        """Deletes an run element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "run", byid, byposition)

        # Delete in Table files
        fileList = []
        exp = _get_all_children(elem, "react")
        for node in exp:
            partit = _get_first_child(node, "partitions")
            if partit is not None:
                finalFileName = "partitions/" + _get_first_child_text(partit, "endPtTable")
                if finalFileName != "partitions/":
                    fileList.append(finalFileName)
        if len(fileList) > 0:
            if self._rdmlFilename is not None and self._rdmlFilename != "":
                if zipfile.is_zipfile(self._rdmlFilename):
                    with zipfile.ZipFile(self._rdmlFilename, 'r') as RDMLin:
                        tempFolder, tempName = tempfile.mkstemp(dir=os.path.dirname(self._rdmlFilename))
                        os.close(tempFolder)
                        with zipfile.ZipFile(tempName, mode='w', compression=zipfile.ZIP_DEFLATED) as RDMLout:
                            RDMLout.comment = RDMLin.comment
                            for item in RDMLin.infolist():
                                if item.filename not in fileList:
                                    RDMLout.writestr(item, RDMLin.read(item.filename))
                    os.remove(self._rdmlFilename)
                    os.rename(tempName, self._rdmlFilename)

        # Delete the node
        self._node.remove(elem)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        allRuns = self.runs()
        runs = []
        for exp in allRuns:
            runs.append(exp.tojson())
        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["runs"] = runs
        return data


class Run:
    """RDML-Python library

    The run element used to read and edit one run.

    Attributes:
        _node: The run node of the RDML XML object.
        _rdmlFilename: The RDML filename.
    """

    def __init__(self, node, rdmlFilename):
        """Inits an run instance.

        Args:
            self: The class self parameter.
            node: The sample node.
            rdmlFilename: The RDML filename.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlFilename = rdmlFilename

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the run subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key in ["description", "instrument", "backgroundDeterminationMethod", "cqDetectionMethod", "runDate"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key == "thermalCyclingConditions":
            forId = _get_first_child(self._node, "thermalCyclingConditions")
            if forId is not None:
                return forId.attrib['id']
            else:
                return None
        if key in ["dataCollectionSoftware_name", "dataCollectionSoftware_version"]:
            ele = _get_first_child(self._node, "dataCollectionSoftware")
            if ele is None:
                return None
            if key == "dataCollectionSoftware_name":
                return _get_first_child_text(ele, "name")
            if key == "dataCollectionSoftware_version":
                return _get_first_child_text(ele, "version")
            raise RdmlError('Run dataCollectionSoftware programming read error.')
        if key in ["pcrFormat_rows", "pcrFormat_columns", "pcrFormat_rowLabel", "pcrFormat_columnLabel"]:
            ele = _get_first_child(self._node, "pcrFormat")
            if ele is None:
                return None
            if key == "pcrFormat_rows":
                return _get_first_child_text(ele, "rows")
            if key == "pcrFormat_columns":
                return _get_first_child_text(ele, "columns")
            if key == "pcrFormat_rowLabel":
                return _get_first_child_text(ele, "rowLabel")
            if key == "pcrFormat_columnLabel":
                return _get_first_child_text(ele, "columnLabel")
            raise RdmlError('Run pcrFormat programming read error.')
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the run subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "cqDetectionMethod":
            if value not in ["", "automated threshold and baseline settings", "manual threshold and baseline settings",
                             "second derivative maximum", "other"]:
                raise RdmlError('Unknown or unsupported run cqDetectionMethod value "' + value + '".')
        if key in ["pcrFormat_rowLabel", "pcrFormat_columnLabel"]:
            if value not in ["ABC", "123", "A1a1"]:
                raise RdmlError('Unknown or unsupported run ' + key + ' value "' + value + '".')

        if key == "id":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["description", "instrument", "backgroundDeterminationMethod", "cqDetectionMethod", "runDate"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key == "thermalCyclingConditions":
            forId = _get_or_create_subelement(self._node, "thermalCyclingConditions", self.xmlkeys())
            if value is not None and value != "":
                # We do not check that ID is valid to allow recreate_lost_ids()
                forId.attrib['id'] = value
            else:
                self._node.remove(forId)
            return
        if key in ["dataCollectionSoftware_name", "dataCollectionSoftware_version"]:
            ele = _get_or_create_subelement(self._node, "dataCollectionSoftware", self.xmlkeys())
            if key == "dataCollectionSoftware_name":
                _change_subelement(ele, "name", ["name", "version"], value, True, "string")
            if key == "dataCollectionSoftware_version":
                _change_subelement(ele, "version", ["name", "version"], value, True, "string")
            _remove_irrelevant_subelement(self._node, "dataCollectionSoftware")
            return
        if key in ["pcrFormat_rows", "pcrFormat_columns", "pcrFormat_rowLabel", "pcrFormat_columnLabel"]:
            ele = _get_or_create_subelement(self._node, "pcrFormat", self.xmlkeys())
            if key == "pcrFormat_rows":
                _change_subelement(ele, "rows", ["rows", "columns", "rowLabel", "columnLabel"], value, True, "string")
            if key == "pcrFormat_columns":
                _change_subelement(ele, "columns", ["rows", "columns", "rowLabel", "columnLabel"], value, True, "string")
            if key == "pcrFormat_rowLabel":
                _change_subelement(ele, "rowLabel", ["rows", "columns", "rowLabel", "columnLabel"], value, True, "string")
            if key == "pcrFormat_columnLabel":
                _change_subelement(ele, "columnLabel", ["rows", "columns", "rowLabel", "columnLabel"], value, True, "string")
            _remove_irrelevant_subelement(self._node, "pcrFormat")
            return
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "instrument", "dataCollectionSoftware_name", "dataCollectionSoftware_version",
                "backgroundDeterminationMethod", "cqDetectionMethod", "thermalCyclingConditions", "pcrFormat_rows",
                "pcrFormat_columns", "pcrFormat_rowLabel", "pcrFormat_columnLabel", "runDate", "react"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "experimenter", "instrument", "dataCollectionSoftware",
                "backgroundDeterminationMethod", "cqDetectionMethod", "thermalCyclingConditions", "pcrFormat",
                "runDate", "react"]

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def experimenter_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "experimenter")

    def update_experimenter_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.experimenter_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "experimenter", id)
                    place = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "experimenter", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_experimenter(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "experimenter", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["experimenters"] = self.experimenter_ids()
        _add_first_child_to_dic(self._node, data, True, "instrument")
        elem = _get_first_child(self._node, "dataCollectionSoftware")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, True, "name")
            _add_first_child_to_dic(elem, qdic, True, "version")
            if len(qdic.keys()) != 0:
                data["dataCollectionSoftware"] = qdic
        _add_first_child_to_dic(self._node, data, True, "backgroundDeterminationMethod")
        _add_first_child_to_dic(self._node, data, True, "cqDetectionMethod")
        forId = _get_first_child(self._node, "thermalCyclingConditions")
        if forId is not None:
            if forId.attrib['id'] != "":
                data["thermalCyclingConditions"] = forId.attrib['id']
        elem = _get_first_child(self._node, "pcrFormat")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "rows")
            _add_first_child_to_dic(elem, qdic, False, "columns")
            _add_first_child_to_dic(elem, qdic, False, "rowLabel")
            _add_first_child_to_dic(elem, qdic, False, "columnLabel")
            data["pcrFormat"] = qdic
        _add_first_child_to_dic(self._node, data, True, "runDate")
        data["react"] = _get_number_of_children(self._node, "react")
        return data

    def export_table(self, dMode):
        """Returns a tab seperated table file with the react fluorescence data.

        Args:
            self: The class self parameter.
            dMode: amp for amplification data, melt for meltcurve data

        Returns:
            A string with the data.
        """

        samTypeLookup = {}
        tarTypeLookup = {}
        tarDyeLookup = {}
        data = ""

        # Get the information for the lookup dictionaries
        pExp = self._node.getparent()
        pRoot = pExp.getparent()
        samples = _get_all_children(pRoot, "sample")
        for sample in samples:
            if sample.attrib['id'] != "":
                samId = sample.attrib['id']
                forType = _get_first_child_text(sample, "type")
                if forType is not "":
                    samTypeLookup[samId] = forType
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType is not "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None:
                    if forId.attrib['id'] != "":
                        tarDyeLookup[tarId] = forId.attrib['id']

        # Now create the header line
        data += "Well\tSample\tSample Type\tTarget\tTarget Type\tDye\t"
        reacts = _get_all_children(self._node, "react")
        if len(reacts) < 1:
            return ""
        react_datas = _get_all_children(reacts[0], "data")
        if len(react_datas) < 1:
            return ""
        headArr = []
        if dMode == "amp":
            adps = _get_all_children(react_datas[0], "adp")
            for adp in adps:
                headArr.append(_get_first_child_text(adp, "cyc"))
            headArr = sorted(headArr, key=int)
        else:
            mdps = _get_all_children(react_datas[0], "mdp")
            for mdp in mdps:
                headArr.append(_get_first_child_text(mdp, "tmp"))
            headArr = sorted(headArr, key=float, reverse=True)
        for hElem in headArr:
            data += hElem + "\t"
        data += '\n'

        # Now create the data lines
        reacts = _get_all_children(self._node, "react")
        wellData = []
        for react in reacts:
            reactId = react.get('id')
            dataSample = reactId + '\t'
            react_sample = "No Sample"
            react_sample_type = "No Sample Type"
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    react_sample = forId.attrib['id']
                    react_sample_type = samTypeLookup[react_sample]
            dataSample += react_sample + '\t' + react_sample_type
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                dataLine = dataSample
                react_target = "No Target"
                react_target_type = "No Target Type"
                react_target_dye = "No Dye"
                forId = _get_first_child(react_data, "tar")
                if forId is not None:
                    if forId.attrib['id'] != "":
                        react_target = forId.attrib['id']
                        react_target_type = tarTypeLookup[react_target]
                        react_target_dye = tarDyeLookup[react_target]
                dataLine += "\t" + react_target + '\t' + react_target_type + '\t' + react_target_dye
                fluorList = []
                if dMode == "amp":
                    adps = _get_all_children(react_data, "adp")
                    for adp in adps:
                        cyc = _get_first_child_text(adp, "cyc")
                        fluor = _get_first_child_text(adp, "fluor")
                        fluorList.append([cyc, fluor])
                    fluorList = sorted(fluorList, key=_sort_list_int)
                else:
                    mdps = _get_all_children(react_data, "mdp")
                    for mdp in mdps:
                        tmp = _get_first_child_text(mdp, "tmp")
                        fluor = _get_first_child_text(mdp, "fluor")
                        fluorList.append([tmp, fluor])
                    fluorList = sorted(fluorList, key=_sort_list_float)
                for hElem in fluorList:
                    dataLine += "\t" + hElem[1]
                dataLine += '\n'
                wellData.append([reactId, dataLine])
        wellData = sorted(wellData, key=_sort_list_int)
        for hElem in wellData:
            data += hElem[1]
        return data

    def import_table(self, rootEl, filename, dMode):
        """Imports data from a tab seperated table file with react fluorescence data.

        Args:
            self: The class self parameter.
            rootEl: The rdml root element.
            filename: The tab file to open.
            dMode: amp for amplification data, melt for meltcurve data.

        Returns:
            A string with the modifications made.
        """

        ret = ""
        with open(filename, "r") as tfile:
            fileContent = tfile.read()

        newlineFix = fileContent.replace("\r\n", "\n")
        tabLines = newlineFix.split("\n")

        head = tabLines[0].split("\t")
        if head[0] != "Well" or head[1] != "Sample" or head[2] != "Sample Type" or head[3] != "Target" \
            or head[4] != "Target Type" or head[5] != "Dye":
            raise RdmlError('The tab-format is not valid, essential columns are missing.')

        # Get the information for the lookup dictionaries
        samTypeLookup = {}
        tarTypeLookup = {}
        dyeLookup = {}
        samples = _get_all_children(rootEl._node, "sample")
        for sample in samples:
            if sample.attrib['id'] != "":
                samId = sample.attrib['id']
                forType = _get_first_child_text(sample, "type")
                if forType is not "":
                    samTypeLookup[samId] = forType
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType is not "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None and forId.attrib['id'] != "":
                    dyeLookup[forId.attrib['id']] = 1

        # Process the lines
        for tabLine in tabLines[1:]:
            sLin = tabLine.split("\t")
            if len(sLin) < 7 or sLin[1] == "" or sLin[2] == "" or sLin[3] == "" or sLin[4] == "" or sLin[5] == "":
                continue
            if sLin[1] not in samTypeLookup:
                rootEl.new_sample(sLin[1], sLin[2])
                samTypeLookup[sLin[1]] = sLin[2]
                ret += "Created sample \"" + sLin[1] + "\" with type \"" + sLin[2] + "\"\n"
            if sLin[3] not in tarTypeLookup:
                if sLin[5] not in dyeLookup:
                    rootEl.new_dye(sLin[5])
                    dyeLookup[sLin[5]] = 1
                    ret += "Created dye \"" + sLin[5] + "\"\n"
                rootEl.new_target(sLin[3], sLin[4])
                elem = rootEl.get_target(byid=sLin[3])
                elem["dyeId"] = sLin[5]
                tarTypeLookup[sLin[3]] = sLin[4]
                ret += "Created " + sLin[3] + " with type \"" + sLin[4] + "\" and dye \"" + sLin[5] + "\"\n"

            react = None
            data = None

            # Get the position number if required
            wellPos = sLin[0]
            if re.search("\D\d+", sLin[0]):
                old_letter = ord(re.sub("\d", "", sLin[0]).upper()) - ord("A")
                old_nr = int(re.sub("\D", "", sLin[0]))
                newId = old_nr + old_letter * int(self["pcrFormat_columns"])
                wellPos = str(newId)
            if re.search("\D\d+\D\d+", sLin[0]):
                old_left = re.sub("\D\d+$", "", sLin[0])
                old_left_letter = ord(re.sub("\d", "", old_left).upper()) - ord("A")
                old_left_nr = int(re.sub("\D", "", old_left)) - 1
                old_right = re.sub("^\D\d+", "", sLin[0])
                old_right_letter = ord(re.sub("\d", "", old_right).upper()) - ord("A")
                old_right_nr = int(re.sub("\D", "", old_right))
                newId = old_left_nr * 8 + old_right_nr + old_left_letter * 768 + old_right_letter * 96
                wellPos = str(newId)

            exp = _get_all_children(self._node, "react")
            for node in exp:
                if wellPos == node.attrib['id']:
                    react = node
                    forId = _get_first_child_text(react, "sample")
                    if forId and forId is not "" and forId.attrib['id'] != sLin[1]:
                        ret += "Missmatch: Well " + wellPos + " (" + sLin[0] + ") has sample \"" + forId.attrib['id'] + \
                               "\" in RDML file and sample \"" + sLin[1] + "\" in tab file.\n"
                    break
            if react is None:
                new_node = ET.Element("react", id=wellPos)
                place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                self._node.insert(place, new_node)
                react = new_node
                new_node = ET.Element("sample", id=sLin[1])
                react.insert(0, new_node)

            exp = _get_all_children(react, "data")
            for node in exp:
                forId = _get_first_child(node, "tar")
                if forId is not None and forId.attrib['id'] == sLin[3]:
                    data = node
                    break
            if data is None:
                new_node = ET.Element("data")
                place = _get_tag_pos(react, "data", ["sample", "data", "partitions"], 9999999)
                react.insert(place, new_node)
                data = new_node
                new_node = ET.Element("tar", id=sLin[3])
                place = _get_tag_pos(data, "tar",
                                     ["tar", "cq", "excl", "adp", "mdp", "endPt", "bgFluor", "quantFluor"],
                                     9999999)
                data.insert(place, new_node)

            if dMode == "amp":
                presentAmp = _get_first_child(data, "adp")
                if presentAmp is not None:
                    ret += "Well " + wellPos + " (" + sLin[0] + ") with sample \"" + sLin[1] + " and target \"" + \
                           sLin[3] + "\" has already amplification data, no data were added.\n"
                else:
                    colCount = 6
                    for col in sLin[6:]:
                        new_node = ET.Element("adp")
                        place = _get_tag_pos(data, "adp",
                                             ["tar", "cq", "excl", "adp", "mdp", "endPt", "bgFluor", "quantFluor"],
                                             9999999)
                        data.insert(place, new_node)
                        new_sub = ET.Element("cyc")
                        new_sub.text = head[colCount]
                        place = _get_tag_pos(new_node, "cyc", ["cyc", "tmp", "fluor"], 9999999)
                        new_node.insert(place, new_sub)
                        new_sub = ET.Element("fluor")
                        new_sub.text = col
                        place = _get_tag_pos(new_node, "fluor", ["cyc", "tmp", "fluor"], 9999999)
                        new_node.insert(place, new_sub)
                        colCount += 1
            if dMode == "melt":
                presentAmp = _get_first_child(data, "mdp")
                if presentAmp is not None:
                    ret += "Well " + wellPos + " (" + sLin[0] + ")  with sample \"" + sLin[1] + " and target \"" + \
                           sLin[3] + "\" has already melting data, no data were added.\n"
                else:
                    colCount = 6
                    for col in sLin[6:]:
                        new_node = ET.Element("mdp")
                        place = _get_tag_pos(data, "mdp",
                                             ["tar", "cq", "excl", "adp", "mdp", "endPt", "bgFluor", "quantFluor"],
                                             9999999)
                        data.insert(place, new_node)
                        new_sub = ET.Element("tmp")
                        new_sub.text = head[colCount]
                        place = _get_tag_pos(new_node, "tmp", ["tmp", "fluor"], 9999999)
                        new_node.insert(place, new_sub)
                        new_sub = ET.Element("fluor")
                        new_sub.text = col
                        place = _get_tag_pos(new_node, "fluor", ["tmp", "fluor"], 9999999)
                        new_node.insert(place, new_sub)
                        colCount += 1
        return ret

    def import_digital_data(self, rootEl, fileformat, filename, filelist):
        """Imports data from a tab seperated table file with digital PCR overview data.

        Args:
            self: The class self parameter.
            rootEl: The rdml root element.
            fileformat: The format of the files (RDML, BioRad).
            filename: The tab overvie file to open (recommended but optional).
            filelist: A list of tab files with fluorescence data (optional, works without filename).

        Returns:
            A string with the modifications made.
        """

        ret = ""
        wellNames = []
        uniqueFileNames = []
        if filelist is None:
            filelist = []

        # Get the information for the lookup dictionaries
        samTypeLookup = {}
        tarTypeLookup = {}
        dyeLookup = {}
        headerLookup = {}
        fileLookup = {}
        fileNameSuggLookup = {}

        samples = _get_all_children(rootEl._node, "sample")
        for sample in samples:
            if sample.attrib['id'] != "":
                samId = sample.attrib['id']
                forType = _get_first_child_text(sample, "type")
                if forType is not "":
                    samTypeLookup[samId] = forType
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType is not "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None and forId.attrib['id'] != "":
                    dyeLookup[forId.attrib['id']] = 1

        # Work the overview file
        if filename is not None:
            with open(filename, "r") as tfile:
                fileContent = tfile.read()

                newlineFix = fileContent.replace("\r\n", "\n")
                tabLines = newlineFix.split("\n")

                posCount = 0

                posWell = 0
                posSample = -1
                posSampleType = -1
                posDye = -1
                posTarget = -1
                posTargetType = -1
                posCopConc = -1
                posPositives = -1
                posNegatives = -1
                posCopConcCh2 = -1
                posPositivesCh2 = -1
                posNegativesCh2 = -1
                posCopConcCh3 = -1
                posPositivesCh3 = -1
                posNegativesCh3 = -1
                posUndefined = -1
                posExcluded = -1
                posVolume = -1
                posFilename = -1

                if fileformat == "RDML":
                    head = tabLines[0].split("\t")
                    for hInfo in head:
                        if hInfo == "Sample":
                            posSample = posCount
                        if hInfo == "SampleType":
                            posSampleType = posCount
                        if hInfo == "Target":
                            posTarget = posCount
                        if hInfo == "TargetType":
                            posTargetType = posCount
                        if hInfo == "Dye":
                            posDye = posCount
                        if hInfo == "Copies":
                            posCopConc = posCount
                        if hInfo == "Positives":
                            posPositives = posCount
                        if hInfo == "Negatives":
                            posNegatives = posCount
                        if hInfo == "Undefined":
                            posUndefined = posCount
                        if hInfo == "Excluded":
                            posExcluded = posCount
                        if hInfo == "Volume":
                            posVolume = posCount
                        if hInfo == "FileName":
                            posFilename = posCount
                        posCount += 1
                elif fileformat == "Bio-Rad":
                    head = tabLines[0].split(";")
                    for hInfo in head:
                        if hInfo == "Sample":
                            posSample = posCount
                        if hInfo in ["TargetType", "TypeAssay"]:
                            posDye = posCount
                        if hInfo in ["Target", "Assay"]:
                            posTarget = posCount
                        if hInfo == "CopiesPer20uLWell":
                            posCopConc = posCount
                        if hInfo == "Positives":
                            posPositives = posCount
                        if hInfo == "Negatives":
                            posNegatives = posCount
                        posCount += 1
                elif fileformat == "Stilla":
                    posWell = 1
                    head = tabLines[0].split(",")
                    for hInfo in head:
                        hInfo = re.sub(r"^ +", '', hInfo)
                        if hInfo == "SampleName":
                            posSample = posCount
                        if hInfo == "Blue_Channel_Concentration":
                            posCopConc = posCount
                        if hInfo == "Blue_Channel_NumberOfPositiveDroplets":
                            posPositives = posCount
                        if hInfo == "Blue_Channel_NumberOfNegativeDroplets":
                            posNegatives = posCount
                        if hInfo == "Green_Channel_Concentration":
                            posCopConcCh2 = posCount
                        if hInfo == "Green_Channel_NumberOfPositiveDroplets":
                            posPositivesCh2 = posCount
                        if hInfo == "Green_Channel_NumberOfNegativeDroplets":
                            posNegativesCh2 = posCount
                        if hInfo == "Red_Channel_Concentration":
                            posCopConcCh3 = posCount
                        if hInfo == "Red_Channel_NumberOfPositiveDroplets":
                            posPositivesCh3 = posCount
                        if hInfo == "Red_Channel_NumberOfNegativeDroplets":
                            posNegativesCh3 = posCount
                        posCount += 1
                    for chan in ["Ch1", "Ch2", "Ch3"]:
                        crTarName = "Target " + chan
                        if crTarName not in tarTypeLookup:
                            if chan not in dyeLookup:
                                rootEl.new_dye(chan)
                                dyeLookup[chan] = 1
                                ret += "Created dye \"" + chan + "\"\n"
                            rootEl.new_target(crTarName, "toi")
                            elem = rootEl.get_target(byid=crTarName)
                            elem["dyeId"] = chan
                            tarTypeLookup[crTarName] = "toi"
                            ret += "Created " + crTarName + " with type \"toi\" and dye \"" + chan + "\"\n"
                else:
                    raise RdmlError('Unknown digital file format.')

                if posSample == -1:
                    raise RdmlError('The overview tab-format is not valid, sample columns are missing.')
                if posDye == -1 and fileformat != "Stilla":
                    raise RdmlError('The overview tab-format is not valid, dye / channel columns are missing.')
                if posTarget == -1 and fileformat != "Stilla":
                    raise RdmlError('The overview tab-format is not valid, target columns are missing.')
                if posPositives == -1:
                    raise RdmlError('The overview tab-format is not valid, positives columns are missing.')
                if posNegatives == -1:
                    raise RdmlError('The overview tab-format is not valid, negatives columns are missing.')

                # Process the lines
                for tabLine in tabLines[1:]:
                    if fileformat == "RDML":
                        sLin = tabLine.split("\t")
                        emptyLine = tabLine.replace("\t", "")
                    elif fileformat == "Bio-Rad":
                        sLin = tabLine.split(";")
                        emptyLine = tabLine.replace(";", "")
                    elif fileformat == "Stilla":
                        tabLine = re.sub(r'^ +', '', tabLine)
                        tabLine = re.sub(r', +', ',', tabLine)
                        sLin = tabLine.split(",")
                        emptyLine1 = tabLine.replace(",", "")
                        emptyLine = emptyLine1.replace(" ", "")
                    else:
                        sLin = tabLine.split("\t")
                        emptyLine = tabLine.replace("\t", "")

                    if len(sLin) < 7 or len(emptyLine) < 5:
                        continue
                    if sLin[posSample] not in samTypeLookup:
                        posSampleTypeName = "unkn"
                        if posSampleType != -1:
                            posSampleTypeName = sLin[posSampleType]
                        rootEl.new_sample(sLin[posSample], posSampleTypeName)
                        samTypeLookup[sLin[posSample]] = posSampleTypeName
                        ret += "Created sample \"" + sLin[posSample] + "\" with type \"" + posSampleTypeName + "\"\n"

                    if fileformat == "RDML":
                        posDyeName = sLin[posDye]
                    elif fileformat == "Bio-Rad":
                        posDyeName = sLin[posDye][:3]
                    elif fileformat == "Stilla":
                        posDyeName = "Not required"
                    else:
                        posDyeName = sLin[posDye]

                    if posTarget != -1:
                        if sLin[posTarget] not in tarTypeLookup:
                            if posDyeName not in dyeLookup:
                                rootEl.new_dye(posDyeName)
                                dyeLookup[posDyeName] = 1
                                ret += "Created dye \"" + posDyeName + "\"\n"
                            posTargetTypeName = "toi"
                            if posTargetType != -1:
                                posTargetTypeName = sLin[posTargetType]
                            rootEl.new_target(sLin[posTarget], posTargetTypeName)
                            elem = rootEl.get_target(byid=sLin[posTarget])
                            elem["dyeId"] = posDyeName
                            tarTypeLookup[sLin[posTarget]] = posTargetTypeName
                            ret += "Created " + sLin[posTarget] + " with type \"" + posTargetTypeName + "\" and dye \"" + posDyeName + "\"\n"

                        if sLin[posWell].upper() not in headerLookup:
                            headerLookup[sLin[posWell].upper()] = {}
                        headerLookup[sLin[posWell].upper()][posDyeName] = sLin[posTarget]

                    if posFilename != -1 and sLin[posFilename] != "":
                        fileNameSuggLookup[sLin[posWell].upper()] = sLin[posFilename]

                    react = None
                    partit = None
                    data = None

                    # Get the position number if required
                    wellPos = sLin[posWell]
                    if fileformat == "Stilla":
                        wellPos = re.sub(r'^\d+-', '', wellPos)

                    if re.search(r"\D\d+", wellPos):
                        old_letter = ord(re.sub(r"\d", "", wellPos.upper())) - ord("A")
                        old_nr = int(re.sub(r"\D", "", wellPos))
                        newId = old_nr + old_letter * int(self["pcrFormat_columns"])
                        wellPos = str(newId)

                    exp = _get_all_children(self._node, "react")
                    for node in exp:
                        if wellPos == node.attrib['id']:
                            react = node
                            forId = _get_first_child_text(react, "sample")
                            if forId and forId is not "" and forId.attrib['id'] != sLin[posSample]:
                                ret += "Missmatch: Well " + wellPos + " (" + sLin[posWell] + ") has sample \"" + forId.attrib['id'] + \
                                       "\" in RDML file and sample \"" + sLin[posSample] + "\" in tab file.\n"
                            break
                    if react is None:
                        new_node = ET.Element("react", id=wellPos)
                        place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                        self._node.insert(place, new_node)
                        react = new_node
                        new_node = ET.Element("sample", id=sLin[posSample])
                        react.insert(0, new_node)

                    partit = _get_first_child(react, "partitions")
                    if partit is None:
                        new_node = ET.Element("partitions")
                        place = _get_tag_pos(react, "partitions", ["sample", "data", "partitions"], 9999999)
                        react.insert(place, new_node)
                        partit = new_node
                        new_node = ET.Element("volume")
                        if fileformat == "RDML":
                            new_node.text = sLin[posVolume]
                        elif fileformat == "Bio-Rad":
                            new_node.text = "0.85"
                        elif fileformat == "Stilla":
                            new_node.text = "0.59"
                        else:
                            new_node.text = "0.70"
                        place = _get_tag_pos(partit, "volume", ["volume", "endPtTable", "data"], 9999999)
                        partit.insert(place, new_node)

                    if fileformat == "Stilla":
                        exp = _get_all_children(partit, "data")
                        for i in range(1, 4):
                            data = None
                            stillaTarget = "Target Ch" + str(i)
                            posDyeName = "Ch" + str(i)
                            stillaConc = "0"
                            stillaPos = "0"
                            stillaNeg = "0"
                            if i == 1:
                                stillaConc = sLin[posCopConc]
                                stillaPos = sLin[posPositives]
                                stillaNeg = sLin[posNegatives]
                            if i == 2:
                                stillaConc = sLin[posCopConcCh2]
                                stillaPos = sLin[posPositivesCh2]
                                stillaNeg = sLin[posNegativesCh2]
                            if i == 3:
                                stillaConc = sLin[posCopConcCh3]
                                stillaPos = sLin[posPositivesCh3]
                                stillaNeg = sLin[posNegativesCh3]

                            if re.search(r"\.", stillaConc):
                                stillaConc = re.sub(r"0+$", "", stillaConc)
                                stillaConc = re.sub(r"\.$", ".0", stillaConc)

                            wellName = re.sub(r'^\d+-', '', sLin[posWell])
                            if wellName.upper() not in headerLookup:
                                headerLookup[wellName.upper()] = {}
                            headerLookup[wellName.upper()][posDyeName] = stillaTarget

                            for node in exp:
                                forId = _get_first_child(node, "tar")
                                if forId is not None and forId.attrib['id'] == stillaTarget:
                                    data = node
                                    break

                            if data is None:
                                new_node = ET.Element("data")
                                place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                                partit.insert(place, new_node)
                                data = new_node
                                new_node = ET.Element("tar", id=stillaTarget)
                                place = _get_tag_pos(data, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                data.insert(place, new_node)

                            new_node = ET.Element("pos")
                            new_node.text = stillaPos
                            place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

                            new_node = ET.Element("neg")
                            new_node.text = stillaNeg
                            place = _get_tag_pos(data, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

                            new_node = ET.Element("conc")
                            new_node.text = stillaConc
                            place = _get_tag_pos(data, "conc", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)
                    else:
                        exp = _get_all_children(partit, "data")
                        for node in exp:
                            forId = _get_first_child(node, "tar")
                            if forId is not None and forId.attrib['id'] == sLin[posTarget]:
                                data = node
                                break
                        if data is None:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            data = new_node
                            new_node = ET.Element("tar", id=sLin[posTarget])
                            place = _get_tag_pos(data, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

                        new_node = ET.Element("pos")
                        new_node.text = sLin[posPositives]
                        place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                        data.insert(place, new_node)

                        new_node = ET.Element("neg")
                        new_node.text = sLin[posNegatives]
                        place = _get_tag_pos(data, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                        data.insert(place, new_node)

                        if posUndefined != -1 and sLin[posUndefined] != "":
                            new_node = ET.Element("undef")
                            new_node.text = sLin[posUndefined]
                            place = _get_tag_pos(data, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

                        if posExcluded != -1 and sLin[posExcluded] != "":
                            new_node = ET.Element("excl")
                            new_node.text = sLin[posExcluded]
                            place = _get_tag_pos(data, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

                        if posCopConc != -1:
                            new_node = ET.Element("conc")
                            if int(sLin[posPositives]) == 0:
                                new_node.text = "0"
                            else:
                                if fileformat == "RDML":
                                    new_node.text = sLin[posCopConc]
                                elif fileformat == "Bio-Rad":
                                    new_node.text = str(float(sLin[posCopConc])/20)
                                else:
                                    new_node.text = sLin[posCopConc]
                            place = _get_tag_pos(data, "conc", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)

        # Read the raw data files

        # Extract the well position from file names
        constNameChars = 0
        if len(filelist) > 0:
            charStopCount = False
            for i in range(len(filelist[0])):
                currChar = None
                if charStopCount is False:
                    for wellFileName in filelist:
                        if currChar is None:
                            currChar = wellFileName[i]
                        else:
                            if currChar != wellFileName[i]:
                                charStopCount = True
                    if charStopCount is False:
                        constNameChars = i + 1

        for wellFileName in filelist:
            currName = wellFileName[constNameChars:].upper()
            currName = currName.replace(".CSV", "")
            currName = currName.replace(".TSV", "")
            currName = currName.replace("_AMPLITUDE", "")
            currName = currName.replace("_COMPENSATEDDATA", "")
            currName = currName.replace("_RAWDATA", "")
            currName = re.sub(r"^\d+_", "", currName)
            wellNames.append(currName)
            fileLookup[currName] = wellFileName

        # Propose a filename for raw data
        runId = self._node.get('id')
        runFix = re.sub(r"[^A-Za-z0-9]", "", runId)
        experimentId = self._node.getparent().get('id')
        experimentFix = re.sub(r"[^A-Za-z0-9]", "", experimentId)
        propFileName = "partitions/" + experimentFix + "_" + runFix

        # Get the used unique file names
        if zipfile.is_zipfile(self._rdmlFilename):
            with zipfile.ZipFile(self._rdmlFilename, 'r') as rdmlObj:
                # Get list of files names in rdml zip
                allRDMLfiles = rdmlObj.namelist()
                for ele in allRDMLfiles:
                    if re.search("^partitions/", ele):
                        uniqueFileNames.append(ele.lower())

        # Now process the files
        warnVolume = ""
        for well in wellNames:
            outTabFile = ""
            keepCh1 = False
            keepCh2 = False
            keepCh3 = False
            header = ""

            react = None
            partit = None
            dataCh1 = None
            dataCh2 = None
            dataCh3 = None

            wellPos = well
            if re.search("\D\d+", well):
                old_letter = ord(re.sub("\d", "", well).upper()) - ord("A")
                old_nr = int(re.sub("\D", "", well))
                newId = old_nr + old_letter * int(self["pcrFormat_columns"])
                wellPos = str(newId)

            exp = _get_all_children(self._node, "react")
            for node in exp:
                if wellPos == node.attrib['id']:
                    react = node
                    break

            if react is None:
                sampleName = "Sample in " + well
                if sampleName not in samTypeLookup:
                    rootEl.new_sample(sampleName, "unkn")
                    samTypeLookup[sampleName] = "unkn"
                    ret += "Created sample \"" + sampleName + "\" with type \"" + "unkn" + "\"\n"
                new_node = ET.Element("react", id=wellPos)
                place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                self._node.insert(place, new_node)
                react = new_node
                new_node = ET.Element("sample", id=sampleName)
                react.insert(0, new_node)

            partit = _get_first_child(react, "partitions")
            if partit is None:
                new_node = ET.Element("partitions")
                place = _get_tag_pos(react, "partitions", ["sample", "data", "partitions"], 9999999)
                react.insert(place, new_node)
                partit = new_node
                new_node = ET.Element("volume")
                if fileformat == "RDML":
                    new_node.text = "0.7"
                    warnVolume = "No information on partition volume given, used 0.7."
                elif fileformat == "Bio-Rad":
                    new_node.text = "0.85"
                elif fileformat == "Stilla":
                    new_node.text = "0.59"
                else:
                    new_node.text = "0.85"
                place = _get_tag_pos(partit, "volume", ["volume", "endPtTable", "data"], 9999999)
                partit.insert(place, new_node)

            if wellPos in fileNameSuggLookup:
                finalFileName = "partitions/" + fileNameSuggLookup[wellPos]
            else:
                finalFileName = "partitions/" + _get_first_child_text(partit, "endPtTable")
                if finalFileName == "partitions/":
                    finalFileName = propFileName + "_" + wellPos + "_" + well + ".tsv"
                    triesCount = 0
                    if finalFileName.lower() in uniqueFileNames:
                        while triesCount < 100:
                            finalFileName = propFileName + "_" + wellPos + "_" + well + "_" + str(triesCount) + ".tsv"
                            if finalFileName.lower() not in uniqueFileNames:
                                uniqueFileNames.append(finalFileName.lower())
                                break

            # print(finalFileName)

            with open(fileLookup[well], "r") as wellfile:
                wellFileContent = wellfile.read()
                newlineFix = wellFileContent.replace("\r\n", "\n")
                wellLines = newlineFix.split("\n")

                if fileformat == "RDML":
                    _writeFileInRDML(self._rdmlFilename, finalFileName, newlineFix)

                    delElem = _get_first_child(partit, "endPtTable")
                    if delElem is not None:
                        partit.remove(delElem)
                    new_node = ET.Element("endPtTable")
                    new_node.text = re.sub(r'^partitions/', '', finalFileName)
                    place = _get_tag_pos(partit, "endPtTable", ["volume", "endPtTable", "data"], 9999999)
                    partit.insert(place, new_node)

                    header = wellLines[0].split("\t")

                    for col in range(0, len(header), 2):
                        cPos = 0
                        cNeg = 0
                        cUndef = 0
                        cExcl = 0

                        if header[col] != "":
                            targetName = header[col]
                            if targetName not in tarTypeLookup:
                                dye = "Ch" + str(int((col + 1) / 2))
                                if dye not in dyeLookup:
                                    rootEl.new_dye(dye)
                                    dyeLookup[dye] = 1
                                    ret += "Created dye \"" + dye + "\"\n"
                                rootEl.new_target(targetName, "toi")
                                elem = rootEl.get_target(byid=targetName)
                                elem["dyeId"] = dye
                                tarTypeLookup[targetName] = "toi"
                                ret += "Created target " + targetName + " with type \"" + "toi" + "\" and dye \"" + dye + "\"\n"

                            for line in wellLines[1:]:
                                splitLine = line.split("\t")
                                if len(splitLine) - 1 < col + 1:
                                    continue
                                if splitLine[col + 1] == "p":
                                    cPos += 1
                                if splitLine[col + 1] == "n":
                                    cNeg += 1
                                if splitLine[col + 1] == "u":
                                    cUndef += 1
                                if splitLine[col + 1] == "e":
                                    cExcl += 1

                            data = None
                            exp = _get_all_children(partit, "data")
                            for node in exp:
                                forId = _get_first_child(node, "tar")
                                if forId is not None and forId.attrib['id'] == targetName:
                                    data = node
                            if data is None:
                                new_node = ET.Element("data")
                                place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                                partit.insert(place, new_node)
                                data = new_node
                                new_node = ET.Element("tar", id=targetName)
                                place = _get_tag_pos(data, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                data.insert(place, new_node)
                            delElem = _get_first_child(partit, "pos")
                            if delElem is not None:
                                data.remove(delElem)
                            new_node = ET.Element("pos")
                            new_node.text = str(cPos)
                            place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)
                            delElem = _get_first_child(partit, "neg")
                            if delElem is not None:
                                data.remove(delElem)
                            new_node = ET.Element("neg")
                            new_node.text = str(cNeg)
                            place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            data.insert(place, new_node)
                            delElem = _get_first_child(partit, "undef")
                            if delElem is not None:
                                data.remove(delElem)
                            if cExcl > 0:
                                new_node = ET.Element("undef")
                                new_node.text = str(cUndef)
                                place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                data.insert(place, new_node)
                            delElem = _get_first_child(partit, "excl")
                            if delElem is not None:
                                data.remove(delElem)
                            if cExcl > 0:
                                new_node = ET.Element("excl")
                                new_node.text = str(cExcl)
                                place = _get_tag_pos(data, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                data.insert(place, new_node)

                elif fileformat == "Bio-Rad":
                    ch1Pos = "0"
                    ch1Neg = "0"
                    ch1sum = 0
                    ch2Pos = "0"
                    ch2Neg = "0"
                    ch2sum = 0

                    if well in headerLookup:
                        if "Ch1" in headerLookup[well]:
                            keepCh1 = True
                            header += headerLookup[well]["Ch1"] + "\t" + headerLookup[well]["Ch1"] + "\t"
                        if "Ch2" in headerLookup[well]:
                            keepCh2 = True
                            header += headerLookup[well]["Ch2"] + "\t" + headerLookup[well]["Ch2"] + "\t"
                        outTabFile += re.sub(r'\t$', '\n', header)
                    else:
                        headerLookup[well] = {}
                        dyes = ["Ch1", "Ch2"]
                        if len(wellLines) > 1:
                            isThereData = wellLines[1].split(",")
                            ch1Pos = ""
                            ch1Neg = ""
                            ch2Pos = ""
                            ch2Neg = ""
                            if re.search(r"\d", isThereData[0]):
                                keepCh1 = True
                            if len(isThereData) > 1 and re.search(r"\d", isThereData[1]):
                                keepCh2 = True
                            for dye in dyes:
                                if dye not in dyeLookup:
                                    rootEl.new_dye(dye)
                                    dyeLookup[dye] = 1
                                    ret += "Created dye \"" + dye + "\"\n"

                            dyeCount = 0
                            for dye in dyes:
                                dyeCount += 1
                                targetName = "Target in " + well + " " + dye
                                if targetName not in tarTypeLookup:
                                    rootEl.new_target(targetName, "toi")
                                    elem = rootEl.get_target(byid=targetName)
                                    elem["dyeId"] = dye
                                    tarTypeLookup[targetName] = "toi"
                                    ret += "Created target " + targetName + " with type \"" + "toi" + "\" and dye \"" + dye + "\"\n"
                                    headerLookup[well][dye] = targetName
                                if (dyeCount == 1 and keepCh1) or (dyeCount == 2 and keepCh2):
                                    header += targetName + "\t" + targetName + "\t"
                            outTabFile += re.sub(r'\t$', '\n', header)

                    if keepCh1 or keepCh2:
                        exp = _get_all_children(partit, "data")
                        for node in exp:
                            forId = _get_first_child(node, "tar")
                            if keepCh1 and forId is not None and forId.attrib['id'] == headerLookup[well]["Ch1"]:
                                dataCh1 = node
                                ch1Pos = _get_first_child_text(dataCh1, "pos")
                                ch1Neg = _get_first_child_text(dataCh1, "neg")
                                ch1sum += int(ch1Pos) + int(ch1Neg)
                            if keepCh2 and forId is not None and forId.attrib['id'] == headerLookup[well]["Ch2"]:
                                dataCh2 = node
                                ch2Pos = _get_first_child_text(dataCh2, "pos")
                                ch2Neg = _get_first_child_text(dataCh2, "neg")
                                ch2sum += int(ch2Pos) + int(ch2Neg)
                        if dataCh1 is None and keepCh1:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh1 = new_node
                            new_node = ET.Element("tar", id=headerLookup[well]["Ch1"])
                            place = _get_tag_pos(dataCh1, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            dataCh1.insert(place, new_node)
                            ch1Pos = ""
                            ch1Neg = ""
                            ch1sum = 2
                        if dataCh2 is None and keepCh2:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh2 = new_node
                            new_node = ET.Element("tar", id=headerLookup[well]["Ch2"])
                            place = _get_tag_pos(dataCh2, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                            dataCh2.insert(place, new_node)
                            ch2Pos = ""
                            ch2Neg = ""
                            ch2sum = 2
                        if dataCh1 is None and dataCh2 is None:
                            continue
                        if ch1sum < 1 and ch2sum < 1:
                            continue

                        if ch1Pos == "" and ch1Neg == "" and ch2Pos == "" and ch2Neg == "":
                            countPart = 0
                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine[0]) < 2:
                                    continue
                                if keepCh1:
                                    outTabFile += splitLine[0] + "\t" + "u"
                                if keepCh2:
                                    if keepCh1:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[1] + "\t" + "u\n"
                                else:
                                    outTabFile += "\n"
                                countPart += 1
                            if keepCh1:
                                new_node = ET.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = ET.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = ET.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh1, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"], 9999999)
                                dataCh1.insert(place, new_node)
                            if keepCh2:
                                new_node = ET.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = ET.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = ET.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh2, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)
                        else:
                            ch1Arr = []
                            ch2Arr = []
                            ch1Cut = 0
                            ch2Cut = 0
                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine) < 2:
                                    continue
                                if keepCh1:
                                    ch1Arr.append(float(splitLine[0]))
                                if keepCh2:
                                    ch2Arr.append(float(splitLine[1]))

                            if keepCh1:
                                ch1Arr.sort()
                                if 0 < int(ch1Neg) <= len(ch1Arr):
                                    ch1Cut = ch1Arr[int(ch1Neg) - 1]
                            if keepCh2:
                                ch2Arr.sort()
                                if 0 < int(ch2Neg) <= len(ch2Arr):
                                    ch2Cut = ch2Arr[int(ch2Neg) - 1]

                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine) < 2:
                                    continue
                                if keepCh1:
                                    outTabFile += splitLine[0] + "\t"
                                    if float(splitLine[0]) > ch1Cut:
                                        outTabFile += "p"
                                    else:
                                        outTabFile += "n"
                                if keepCh2:
                                    if keepCh1:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[1] + "\t"
                                    if float(splitLine[1]) > ch2Cut:
                                        outTabFile += "p\n"
                                    else:
                                        outTabFile += "n\n"
                                else:
                                    outTabFile += "\n"
                        _writeFileInRDML(self._rdmlFilename, finalFileName, outTabFile)
                        new_node = ET.Element("endPtTable")
                        new_node.text = re.sub(r'^partitions/', '', finalFileName)
                        place = _get_tag_pos(partit, "endPtTable", ["volume", "endPtTable", "data"], 9999999)
                        partit.insert(place, new_node)
                    else:
                        react.remove(partit)
                elif fileformat == "Stilla":
                    ch1Pos = "0"
                    ch1Neg = "0"
                    ch1sum = 0
                    ch2Pos = "0"
                    ch2Neg = "0"
                    ch2sum = 0
                    ch3Pos = "0"
                    ch3Neg = "0"
                    ch3sum = 0

                    if well in headerLookup:
                        if "Ch1" in headerLookup[well]:
                            keepCh1 = True
                            header += headerLookup[well]["Ch1"] + "\t" + headerLookup[well]["Ch1"] + "\t"
                        if "Ch2" in headerLookup[well]:
                            keepCh2 = True
                            header += headerLookup[well]["Ch2"] + "\t" + headerLookup[well]["Ch2"] + "\t"
                        if "Ch3" in headerLookup[well]:
                            keepCh3 = True
                            header += headerLookup[well]["Ch3"] + "\t" + headerLookup[well]["Ch3"] + "\t"
                        outTabFile += re.sub(r'\t$', '\n', header)
                    else:
                        headerLookup[well] = {}
                        dyes = ["Ch1", "Ch2", "Ch3"]
                        if len(wellLines) > 1:
                            isThereData = wellLines[1].split(",")
                            ch1Pos = ""
                            ch1Neg = ""
                            ch2Pos = ""
                            ch2Neg = ""
                            ch3Pos = ""
                            ch3Neg = ""
                            if re.search(r"\d", isThereData[0]):
                                keepCh1 = True
                            if len(isThereData) > 1 and re.search(r"\d", isThereData[1]):
                                keepCh2 = True
                            if len(isThereData) > 2 and re.search(r"\d", isThereData[2]):
                                keepCh3 = True
                            for dye in dyes:
                                if dye not in dyeLookup:
                                    rootEl.new_dye(dye)
                                    dyeLookup[dye] = 1
                                    ret += "Created dye \"" + dye + "\"\n"
                            dyeCount = 0
                            for dye in dyes:
                                dyeCount += 1
                                targetName = "Target in " + well + " " + dye
                                if targetName not in tarTypeLookup:
                                    rootEl.new_target(targetName, "toi")
                                    elem = rootEl.get_target(byid=targetName)
                                    elem["dyeId"] = dye
                                    tarTypeLookup[targetName] = "toi"
                                    ret += "Created target " + targetName + " with type \"" + "toi" + "\" and dye \"" + dye + "\"\n"
                                    if (dyeCount == 1 and keepCh1) or (dyeCount == 2 and keepCh2) or (dyeCount == 3 and keepCh3):
                                        headerLookup[well][dye] = targetName
                                header += targetName + "\t" + targetName + "\t"
                            outTabFile += re.sub(r'\t$', '\n', header)

                    if keepCh1 or keepCh2 or keepCh3:
                        exp = _get_all_children(partit, "data")
                        for node in exp:
                            forId = _get_first_child(node, "tar")
                            if keepCh1 and forId is not None and forId.attrib['id'] == headerLookup[well]["Ch1"]:
                                dataCh1 = node
                                ch1Pos = _get_first_child_text(dataCh1, "pos")
                                ch1Neg = _get_first_child_text(dataCh1, "neg")
                                ch1sum += int(ch1Pos) + int(ch1Neg)
                            if keepCh2 and forId is not None and forId.attrib['id'] == headerLookup[well]["Ch2"]:
                                dataCh2 = node
                                ch2Pos = _get_first_child_text(dataCh2, "pos")
                                ch2Neg = _get_first_child_text(dataCh2, "neg")
                                ch2sum += int(ch2Pos) + int(ch2Neg)
                            if keepCh3 and forId is not None and forId.attrib['id'] == headerLookup[well]["Ch3"]:
                                dataCh3 = node
                                ch3Pos = _get_first_child_text(dataCh3, "pos")
                                ch3Neg = _get_first_child_text(dataCh3, "neg")
                                ch3sum += int(ch3Pos) + int(ch3Neg)
                        if dataCh1 is None and keepCh1:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh1 = new_node
                            new_node = ET.Element("tar", id=headerLookup[well]["Ch1"])
                            place = _get_tag_pos(dataCh1, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                 9999999)
                            dataCh1.insert(place, new_node)
                            ch1Pos = ""
                            ch1Neg = ""
                            ch1sum = 2
                        if dataCh2 is None and keepCh2:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh2 = new_node
                            new_node = ET.Element("tar", id=headerLookup[well]["Ch2"])
                            place = _get_tag_pos(dataCh2, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                 9999999)
                            dataCh2.insert(place, new_node)
                            ch2Pos = ""
                            ch2Neg = ""
                            ch2sum = 2
                        if dataCh3 is None and keepCh3:
                            new_node = ET.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh3 = new_node
                            new_node = ET.Element("tar", id=headerLookup[well]["Ch3"])
                            place = _get_tag_pos(dataCh3, "tar", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                 9999999)
                            dataCh3.insert(place, new_node)
                            ch3Pos = ""
                            ch3Neg = ""
                            ch3sum = 2
                        if dataCh1 is None and dataCh2 is None and dataCh3 is None:
                            continue
                        if ch1sum < 1 and ch2sum < 1 and ch3sum < 1:
                            continue

                        if ch1Pos == "" and ch1Neg == "" and ch2Pos == "" and ch2Neg == "" and ch3Pos == "" and ch3Neg == "":
                            countPart = 0
                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine[0]) < 2:
                                    continue
                                if keepCh1:
                                    outTabFile += splitLine[0] + "\t" + "u"
                                if keepCh2:
                                    if keepCh1:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[1] + "\t" + "u"
                                if keepCh3:
                                    if keepCh1 or keepCh2:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[2] + "\t" + "u\n"
                                else:
                                    outTabFile += "\n"
                                countPart += 1
                            if keepCh1:
                                new_node = ET.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh1.insert(place, new_node)

                                new_node = ET.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh1.insert(place, new_node)

                                new_node = ET.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh1, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh1.insert(place, new_node)
                            if keepCh2:
                                new_node = ET.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = ET.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = ET.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh2, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh2.insert(place, new_node)
                            if keepCh3:
                                new_node = ET.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh3, "pos", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh3.insert(place, new_node)

                                new_node = ET.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh3, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh3.insert(place, new_node)

                                new_node = ET.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh3, "neg", ["tar", "pos", "neg", "undef", "excl", "conc"],
                                                     9999999)
                                dataCh3.insert(place, new_node)
                        else:
                            ch1Arr = []
                            ch2Arr = []
                            ch3Arr = []
                            ch1Cut = 0
                            ch2Cut = 0
                            ch3Cut = 0
                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine) < 3:
                                    continue
                                if keepCh1:
                                    ch1Arr.append(float(splitLine[0]))
                                if keepCh2:
                                    ch2Arr.append(float(splitLine[1]))
                                if keepCh3:
                                    ch3Arr.append(float(splitLine[2]))

                            if keepCh1:
                                ch1Arr.sort()
                                if 0 < int(ch1Neg) <= len(ch1Arr):
                                    ch1Cut = ch1Arr[int(ch1Neg) - 1]
                            if keepCh2:
                                ch2Arr.sort()
                                if 0 < int(ch2Neg) <= len(ch2Arr):
                                    ch2Cut = ch2Arr[int(ch2Neg) - 1]
                            if keepCh3:
                                ch3Arr.sort()
                                if 0 < int(ch3Neg) <= len(ch3Arr):
                                    ch3Cut = ch3Arr[int(ch3Neg) - 1]

                            for line in wellLines[1:]:
                                splitLine = line.split(",")
                                if len(splitLine) < 2:
                                    continue
                                if keepCh1:
                                    outTabFile += splitLine[0] + "\t"
                                    if float(splitLine[0]) > ch1Cut:
                                        outTabFile += "p"
                                    else:
                                        outTabFile += "n"
                                if keepCh2:
                                    if keepCh1:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[1] + "\t"
                                    if float(splitLine[1]) > ch2Cut:
                                        outTabFile += "p"
                                    else:
                                        outTabFile += "n"
                                if keepCh3:
                                    if keepCh1 or keepCh2:
                                        outTabFile += "\t"
                                    outTabFile += splitLine[2] + "\t"
                                    if float(splitLine[2]) > ch3Cut:
                                        outTabFile += "p\n"
                                    else:
                                        outTabFile += "n\n"
                                else:
                                    outTabFile += "\n"
                        _writeFileInRDML(self._rdmlFilename, finalFileName, outTabFile)
                        new_node = ET.Element("endPtTable")
                        new_node.text = re.sub(r'^partitions/', '', finalFileName)
                        place = _get_tag_pos(partit, "endPtTable", ["volume", "endPtTable", "data"], 9999999)
                        partit.insert(place, new_node)
                    else:
                        react.remove(partit)

        ret += warnVolume
        return ret

    def get_digital_overview_data(self, rootEl):
        """Provides the digital overview data in tab seperated format.

        Args:
            self: The class self parameter.
            rootEl: The rdml root element.

        Returns:
            A string with the overview data table.
        """

        #       0    1      2        3          4          5      6      7         8         9         10        11        12        13
        ret = "Pos\tWell\tSample\tSampleType\tTarget\tTargetType\tDye\tCopies\tPositives\tNegatives\tUndefined\tExcluded\tVolume\tFileName\n"
        tabLines = []

        # Fill the lookup dics
        samTypeLookup = {}
        tarTypeLookup = {}
        tarDyeLookup = {}

        samples = _get_all_children(rootEl._node, "sample")
        for sample in samples:
            if sample.attrib['id'] != "":
                samId = sample.attrib['id']
                forType = _get_first_child_text(sample, "type")
                if forType is not "":
                    samTypeLookup[samId] = forType
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType is not "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None and forId.attrib['id'] != "":
                    tarDyeLookup[tarId] = forId.attrib['id']

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            pPos = react.attrib['id']
            posId = int(react.attrib['id'])
            pIdNumber = posId % int(self["pcrFormat_columns"])
            pIdLetter = chr(ord("A") + int(posId / int(self["pcrFormat_columns"])))
            pWell = pIdLetter + str(pIdNumber)
            pSample = ""
            pSampleType = ""
            pFileName = ""
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    pSample = forId.attrib['id']
                    pSampleType = samTypeLookup[forId.attrib['id']]
            partit = _get_first_child(react, "partitions")
            if partit is not None:
                endPtTable = _get_first_child_text(partit, "endPtTable")
                if endPtTable is not "":
                    pFileName = endPtTable
                pVolume = _get_first_child_text(partit, "volume")
                partit_datas = _get_all_children(partit, "data")
                for partit_data in partit_datas:
                    pTarget = ""
                    pTargetType = ""
                    pDye = ""
                    forId = _get_first_child(partit_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] != "":
                            pTarget = forId.attrib['id']
                            pTargetType = tarTypeLookup[pTarget]
                            pDye = tarDyeLookup[pTarget]
                    pCopies = _get_first_child_text(partit_data, "conc")
                    pPositives = _get_first_child_text(partit_data, "pos")
                    pNegatives = _get_first_child_text(partit_data, "neg")
                    pUnknown = _get_first_child_text(partit_data, "undef")
                    pExcluded = _get_first_child_text(partit_data, "excl")

                    retLine = pPos + "\t"
                    retLine += pWell + "\t"
                    retLine += pSample + "\t"
                    retLine += pSampleType + "\t"
                    retLine += pTarget + "\t"
                    retLine += pTargetType + "\t"
                    retLine += pDye + "\t"
                    retLine += pCopies + "\t"
                    retLine += pPositives + "\t"
                    retLine += pNegatives + "\t"
                    retLine += pUnknown + "\t"
                    retLine += pExcluded + "\t"
                    retLine += pVolume + "\t"
                    retLine += pFileName + "\n"
                    tabLines.append(retLine)
        tabLines.sort(key=_sort_list_digital_PCR)
        for tLine in tabLines:
            ret += tLine
        return ret

    def get_digital_raw_data(self, reactPos):
        """Provides the digital of a react in tab seperated format.

        Args:
            self: The class self parameter.
            reactPos: The react id to get the digital raw data from

        Returns:
            A string with the raw data table.
        """

        react = None
        partit = None
        retVal = ""

        # Get the position number if required
        wellPos = str(reactPos)
        if re.search("\D\d+", wellPos):
            old_letter = ord(re.sub("\d", "", wellPos.upper())) - ord("A")
            old_nr = int(re.sub("\D", "", wellPos))
            newId = old_nr + old_letter * int(self["pcrFormat_columns"])
            wellPos = str(newId)

        exp = _get_all_children(self._node, "react")
        for node in exp:
            if wellPos == node.attrib['id']:
                react = node
                break
        if react is None:
            return ""

        partit = _get_first_child(react, "partitions")
        if partit is None:
            return ""

        finalFileName = "partitions/" + _get_first_child_text(partit, "endPtTable")
        if finalFileName == "partitions/":
            return ""

        if zipfile.is_zipfile(self._rdmlFilename):
            zf = zipfile.ZipFile(self._rdmlFilename, 'r')
            try:
                retVal = zf.read(finalFileName).decode('utf-8')
            except KeyError:
                raise RdmlError('No ' + finalFileName + ' in compressed RDML file found.')
            finally:
                zf.close()
        return retVal

    def getreactjson(self):
        """Returns a json of the react data including fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        all_data = {}
        data = []
        reacts = _get_all_children(self._node, "react")

        adp_cyc_max = 0.0
        adp_fluor_min = 99999999
        adp_fluor_max = 0.0
        mdp_tmp_min = 120.0
        mdp_tmp_max = 0.0
        mdp_fluor_min = 99999999
        mdp_fluor_max = 0.0
        max_data = 0
        max_partition_data = 0
        for react in reacts:
            react_json = {
                "id": react.get('id'),
            }
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    react_json["sample"] = forId.attrib['id']
            react_datas = _get_all_children(react, "data")
            max_data = max(max_data, len(react_datas))
            react_datas_json = []
            for react_data in react_datas:
                in_react = {}
                forId = _get_first_child(react_data, "tar")
                if forId is not None:
                    if forId.attrib['id'] != "":
                        in_react["tar"] = forId.attrib['id']
                _add_first_child_to_dic(react_data, in_react, True, "cq")
                _add_first_child_to_dic(react_data, in_react, True, "excl")
                _add_first_child_to_dic(react_data, in_react, True, "endPt")
                _add_first_child_to_dic(react_data, in_react, True, "bgFluor")
                _add_first_child_to_dic(react_data, in_react, True, "bgFluorSlp")
                _add_first_child_to_dic(react_data, in_react, True, "quantFluor")
                adps = _get_all_children(react_data, "adp")
                adps_json = []
                for adp in adps:
                    cyc = _get_first_child_text(adp, "cyc")
                    fluor = _get_first_child_text(adp, "fluor")
                    adp_cyc_max = max(adp_cyc_max, float(cyc))
                    adp_fluor_min = min(adp_fluor_min, float(fluor))
                    adp_fluor_max = max(adp_fluor_max, float(fluor))
                    in_adp = [cyc, fluor, _get_first_child_text(adp, "tmp")]
                    adps_json.append(in_adp)
                in_react["adps"] = adps_json
                mdps = _get_all_children(react_data, "mdp")
                mdps_json = []
                for mdp in mdps:
                    tmp = _get_first_child_text(mdp, "tmp")
                    fluor = _get_first_child_text(mdp, "fluor")
                    mdp_tmp_min = min(mdp_tmp_min, float(tmp))
                    mdp_tmp_max = max(mdp_tmp_max, float(tmp))
                    mdp_fluor_min = min(mdp_fluor_min, float(fluor))
                    mdp_fluor_max = max(mdp_fluor_max, float(fluor))
                    in_mdp = [tmp, fluor]
                    mdps_json.append(in_mdp)
                in_react["mdps"] = mdps_json
                react_datas_json.append(in_react)
            react_json["datas"] = react_datas_json
            partit = _get_first_child(react, "partitions")
            if partit is not None:
                in_partitions = {}
                endPtTable = _get_first_child_text(partit, "endPtTable")
                if endPtTable is not "":
                    in_partitions["endPtTable"] = endPtTable
                partit_datas = _get_all_children(partit, "data")
                max_partition_data = max(max_partition_data, len(partit_datas))
                partit_datas_json = []
                for partit_data in partit_datas:
                    in_partit = {}
                    forId = _get_first_child(partit_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] != "":
                            in_partit["tar"] = forId.attrib['id']
                    _add_first_child_to_dic(partit_data, in_partit, False, "pos")
                    _add_first_child_to_dic(partit_data, in_partit, False, "neg")
                    _add_first_child_to_dic(partit_data, in_partit, True, "undef")
                    _add_first_child_to_dic(partit_data, in_partit, True, "excl")
                    _add_first_child_to_dic(partit_data, in_partit, True, "conc")
                    partit_datas_json.append(in_partit)
                in_partitions["datas"] = partit_datas_json
                react_json["partitions"] = in_partitions
            data.append(react_json)
        all_data["reacts"] = data
        all_data["adp_cyc_max"] = adp_cyc_max
        all_data["adp_fluor_min"] = adp_fluor_min
        all_data["adp_fluor_max"] = adp_fluor_max
        all_data["mdp_tmp_min"] = mdp_tmp_min
        all_data["mdp_tmp_max"] = mdp_tmp_max
        all_data["mdp_fluor_min"] = mdp_fluor_min
        all_data["mdp_fluor_max"] = mdp_fluor_max
        all_data["max_data_len"] = max_data
        all_data["max_partition_data_len"] = max_partition_data
        return all_data

    def linRegPCR(self, perTarget=True, baselineCorr=True, commaConv=False, ignoreExclusion=False,
                  saveRaw=False, saveBaslineCorr=False, saveResultsList=False, saveResultsCSV=False, verbose=False):
        """Performs LinRegPCR on the run. Mofifies the cq values and returns a json with additional data.

        Args:
            self: The class self parameter.
            perTarget: If true, calculate efficiency per target, if false, for all samples.
            baselineCorr: If true, do baseline correction for all samples.
            commaConv: If true, convert comma separator to dot.
            ignoreExclusion: If true, ignore the RDML exclusion strings.
            saveRaw: If true, no raw values are given in the returned data
            saveBaslineCorr: If true, no baseline corrected values are given in the returned data
            saveResultsList: If true, return a 2d array object.
            saveResultsCSV: If true, return a csv string.
            verbose: If true, comment every performed step.

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            rawData: A 2d array with the raw fluorescence values
            baselineCorrectedData: A 2d array with the baseline corrected raw fluorescence values
            resultsList: A 2d array object.
            resultsCSV: A csv string.
        """

        ##############################
        # Collect the data in arrays #
        ##############################

        # res is a 2 dimensional array accessed only by
        # variables, so columns might be added here
        header = [["id",  # 0
                   "sample",  # 1
                   "target",   # 2
                   "excluded",   # 3
                   "baseline",   # 4
                   "lower limit",   # 5
                   "threshold",   # 6
                   "upper limit",   # 7
                   "n included",   # 8
                   "Cq",   # 9
                   "indiv PCR eff",   # 10
                   "N0 (indiv eff)",   # 11
                   "R2",   # 12
                   "mean PCR eff",   # 13
                   "N0 (mean eff)",   # 14
                   "no amplification",   # 15
                   "baseline error",   # 16
                   "no plateau",   # 17
                   "noisy sample",   # 18
                   "PCR efficiency outside 5%"]]   # 19
        rar_id = 0
        rar_sample = 1
        rar_tar = 2
        rar_excl = 3
        rar_baseline = 4
        rar_lower_limit = 5
        rar_threshold = 6
        rar_upper_limit = 7
        rar_n_included = 8
        rar_Cq = 9
        rar_indiv_PCR_eff = 10
        rar_N0_indiv_eff = 11
        rar_R2 = 12
        rar_mean_PCR_eff = 13
        rar_N0_mean_eff = 14
        rar_no_amplification = 15
        rar_baseline_error = 16
        rar_no_plateau = 17
        rar_noisy_sample = 18
        rar_PCR_efficiency_outside = 19

        res = []
        finalData = {}
        adp_cyc_max = 0

        reacts = _get_all_children(self._node, "react")

        # First get the max number of cycles and create the numpy array
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                adps = _get_all_children(react_data, "adp")
                for adp in adps:
                    cyc = _get_first_child_text(adp, "cyc")
                    adp_cyc_max = max(adp_cyc_max, float(cyc))
        adp_cyc_max = math.ceil(adp_cyc_max)

        # spFl is the shape for all fluorescence numpy data arrays
        spFl = (len(reacts), adp_cyc_max)
        rawFluor = np.zeros(spFl, dtype=np.float64)  # Fixme: nan
        # Initialization of the vecNoAmplification vector
        vecExcludedByUser = np.zeros(spFl[0], dtype=np.bool)
        vecTarget = np.zeros(spFl[0], dtype=np.int)
        vecTarget[vecTarget <= 0] = -1

        # Now process the data for numpy and create results array
        rowCount = 0
        for react in reacts:
            posId = react.get('id')
            sample = ""
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    sample = forId.attrib['id']
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                forId = _get_first_child(react_data, "tar")
                target = ""
                if forId is not None:
                    if forId.attrib['id'] != "":
                        target = forId.attrib['id']
                # _add_first_child_to_dic(react_data, in_react, True, "cq")
                if ignoreExclusion:
                    excl = ""
                else:
                    excl = _get_first_child_text(react_data, "excl")
                if not excl == "":
                    vecExcludedByUser[rowCount] = True
                res.append([posId, sample, target, excl, "", "", "", "", "", "",
                            "", "", "", "", "", "", "", "", "", "", ""])
                adps = _get_all_children(react_data, "adp")
                for adp in adps:
                    cyc = int(math.ceil(float(_get_first_child_text(adp, "cyc")))) - 1
                    fluor = _get_first_child_text(adp, "fluor")
                    if commaConv:
                        noDot = fluor.replace(".", "")
                        fluor = noDot.replace(",", ".")
                    rawFluor[rowCount, cyc] = float(fluor)
                rowCount += 1

        if saveRaw:
            rawArr = [[header[0][rar_id], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl]]]
            for i in range(0, spFl[1]):
                rawArr[0].append(i + 1)
            for i in range(0, spFl[0]):
                rawArr.append([res[i][rar_id], res[i][rar_sample], res[i][rar_tar], res[i][rar_excl]])
                for k in range(0, spFl[1]):
                    rawArr[i + 1].append(float(rawFluor[i, k]))
            finalData["rawData"] = rawArr

        # Create the target variables
        # Position 0 is for the general over all window without targets
        targetsCount = 1
        tarWinLookup = {}
        for i in range(0, spFl[0]):
            if res[i][2] not in tarWinLookup:
                tarWinLookup[res[i][2]] = targetsCount
                vecTarget[i] = targetsCount
                targetsCount += 1
            vecTarget[i] = tarWinLookup[res[i][2]]
        upwin = np.zeros(targetsCount, dtype=np.float64)
        lowwin = np.zeros(targetsCount, dtype=np.float64)

        ########################
        # Baseline correction  #
        ########################
        # Fixme: Delete
        start_time = dt.datetime.now()
        outFFile = ""

        # Initialization of the vecNoAmplification vector
        vecNoAmplification = np.zeros(spFl[0], dtype=np.bool)

        # This is the vector for the baseline error quality check
        vecNoPlateau = np.zeros(spFl[0], dtype=np.bool)

        # This is the vector for the baseline error quality check
        vecBaselineError = np.zeros(spFl[0], dtype=np.bool)

        # This is the vector for the skip sample check
        vecSkipSample = np.zeros(spFl[0], dtype=np.bool)

        vecShortloglin = np.zeros(spFl[0], dtype=np.bool)

        vecIsUsedInWoL = np.zeros(spFl[0], dtype=np.bool)

        pcreff = np.zeros(spFl[0], dtype=np.float64)

        if baselineCorr:
            ####################
            # SDM calculation  #
            ####################

            # First and Second Derivative values calculation

            # Shifted matrix of the raw data
            rawFluorShift = np.roll(rawFluor, 1, axis=1)  # Shift to right - real position is -0.5
            rawFluorShift[:, 0] = np.nan
            # Subtraction of the shifted matrix to the raw data
            firstDerivative = rawFluor - rawFluorShift

            # Shifted matrix of the firstDerivative
            firstDerivativeShift = np.roll(firstDerivative, -1, axis=1)  # Shift to left
            firstDerivativeShift[:, -1] = np.nan
            # Subtraction of the firstDerivative values to the shifted matrix
            secondDerivative = firstDerivativeShift - firstDerivative

            # Calculation of the second derivative maximum per well
            # SDM is the max flo value per well and SDMcycles to the corresponding cycle
            FDM = np.nanmax(firstDerivative, axis=1)
            FDMcycles = np.nanargmax(firstDerivative, axis=1) + 1
            SDM = np.nanmax(secondDerivative, axis=1)
            SDMcycles = np.nanargmax(secondDerivative, axis=1) + 1
            _numpyTwoAxisSave(secondDerivative, "linPas/ot_secondDerivative.tsv")

            ########################################
            # Start point of the exponential phase #
            ########################################

            # For loop which determine the start of the exponential phase

            # This first loop creates the vector that will be used during the
            # initialisation phase of the baseline estimation

            # Initialisation
            startExpPhase = np.zeros(spFl[0], dtype=np.float64)
            startExpCycles = np.zeros(spFl[0], dtype=np.int)

            # for each row of the data
            for i in range(0, spFl[0]):
                # the loop starts from the SDM cycle of each row
                j = SDMcycles[i]
                while rawFluor[i, j] >= rawFluor[i, j - 1] and j > 1:
                    j = j - 1
                startExpPhase[i] = rawFluor[i, j]
                startExpCycles[i] = j + 1

            # For loop which modifies the start cycle if too far from the SDM cycle

            # This loop looks for the cases where the previously calculated start point cycle is
            # more than 10 cycles far from the SDM point. When it finds such a case, it replaces
            # the start cycle by the SDM cycle minus 10. This vector will be used after the initialization
            # phase

            # Initialization of the expoPhaseBaselineEstimation vector
            expoPhaseBaselineEstimation = startExpCycles.copy()

     #       expoPhaseBaselineEstimation[SDMcycles - startExpCycles > 10] = SDMcycles - 10
            for i in range(0, startExpCycles.shape[0]):
                if SDMcycles[i] - startExpCycles[i] > 10:
                    expoPhaseBaselineEstimation[i] = SDMcycles[i] - 10

            ###########################################################################
            # First quality check : Is there enough amplification during the reaction #
            ###########################################################################

            # Slope calculation per well

            logical = np.ones(spFl)
            vecCycle2 = np.arange(1, (spFl[1] + 1), dtype=np.int)

            # the intercept is never used in the code but has to be mentioned otherwise
            # the linearRegression can't be called.

            rawMod = rawFluor.copy()
            rawMod[np.isnan(rawMod)] = 0
            rawMod[rawMod <= 0.00000001] = np.nan
            [slopeAmp, interceptAmp] = _linearRegression(vecCycle2[np.newaxis, :], np.log10(rawMod), logical)

            # Minimum of fluorescence values per well
            minFlu = np.nanmin(rawMod, axis=1)
            sampmin = minFlu.copy()
            bgarr = 0.99 * minFlu
            defbgarr = bgarr.copy()
            cycleMinFlu = np.nanargmin(rawMod, axis=1) + 1
            minFluMat = rawMod - bgarr[:, np.newaxis]

            minFluMat[np.isnan(minFluMat)] = 0
            minFluMat[minFluMat <= 0.00000001] = np.nan
            minFluMatCount = minFluMat.copy()
            minFluMatMean = minFluMat.copy()
            for i in range(0, spFl[0]):
                minFluMatMean[i, cycleMinFlu[i] - 1] = np.nan
            minFluMatCount[np.isnan(minFluMatCount)] = 0
            minFluMatCount[minFluMatCount > 0] = 1
            minFluMatCountSum = np.sum(minFluMatCount, axis=1)

            [minSlopeAmp, interceptAmp] = _linearRegression(vecCycle2[np.newaxis, :], np.log10(minFluMat), logical)

            # Check to detect the negative slopes and the PCR reactions that have an
            # amplification less than seven the minimum fluorescence

            fstop = FDMcycles.copy()  # Fixme: only for correct length
            fstop2 = FDMcycles.copy()  # Fixme: only for correct length
            fstart = FDMcycles.copy()  # Fixme: only for correct length
            fstart2 = FDMcycles.copy()  # Fixme: only for correct length

            for i in range(0, spFl[0]):  # Fixme: use np.nanmean with isnan check
                if slopeAmp[i] < 0 or minSlopeAmp[i] < (np.log10(7.0) / minFluMatCountSum[i]) or \
                   ((((minFluMat[i, 8] + minFluMat[i, 9]) / 2) / ((minFluMat[i, 0] + minFluMat[i, 1]) / 2) < 1.2) and
                    (minFluMat[i, -1] / np.nanmean(minFluMatMean[i, 0:10])) < 7):
                    vecNoAmplification[i] = True

                if not vecNoAmplification[i]:
                    fstop[i] = _findLRstop(minFluMat, i)
                    [fstart[i], fstart2[i]] = _findLRstart(minFluMat, i, fstop[i])
                    fstop2[i] = fstop[i]
                else:
                    vecSkipSample[i] = True
                    fstop[i] = minFluMat.shape[1]
                    fstart[i] = 1
                    fstart2[i] = 1

                if fstop[i] == minFluMat.shape[1]:
                    vecNoPlateau[i] = True

                if minFluMat[i, fstop[i] - 1] <= minFluMat[i, fstop[i] - 2] <= minFluMat[i, fstop[i] - 3]:
                    vecNoAmplification[i] = True
                    vecSkipSample[i] = True

            # Set an initial window
            meanmax = _GetMeanMax(minFluMat, None, vecSkipSample, vecNoPlateau)
            upwin[0] = np.log10(0.1 * meanmax)
            lowwin[0] = np.log10(0.1 * meanmax / 16.0)

            _numpyTwoAxisSave(vecNoAmplification, "linPas/numpy_vecNoAmplification.tsv")
            _numpyTwoAxisSave(vecNoPlateau, "linPas/numpy_vecNoPlateau.tsv")
            _numpyTwoAxisSave(minFluMat, "linPas/numpy_fluor_data_1.tsv")
            _numpyTwoAxisSave(fstop, "linPas/numpy_fstop_1.tsv")
            _numpyTwoAxisSave(fstart, "linPas/numpy_fstart_1.tsv")
          #  for i in range(0, spFl[0]):
           #     fstop[i] = _findLRstop(minFluMat, i)
           #     [ddldldl, fstart[i]] = _findLRstart(minFluMat, i, fstop[i])
            _numpyTwoAxisSave(fstop2, "linPas/numpy_fstop_2.tsv")
            _numpyTwoAxisSave(fstart2, "linPas/numpy_fstart_2.tsv")

            # Consecutive derivative points near the SDM (caution arrays start at 0, not 1!)
 #           for i in range(0, firstDerivativeShift.shape[0]):
  #              if not (firstDerivativeShift[i, SDMcycles[i] - 4] <= firstDerivativeShift[i, SDMcycles[i] - 3] <=
   #                     firstDerivativeShift[i, SDMcycles[i] - 2] <= firstDerivativeShift[i, SDMcycles[i] - 1] or
    #                    firstDerivativeShift[i, SDMcycles[i] - 3] <= firstDerivativeShift[i, SDMcycles[i] - 2] <=
     #                   firstDerivativeShift[i, SDMcycles[i] - 1] <= firstDerivativeShift[i, SDMcycles[i]]):
      #              vecNoAmplification[i] = True

            ##################################################
            # Main loop : Calculation of the baseline values #
            ##################################################

            # The output vector of the main for loops need to be initialized before the start of the loop :
            # This will be the baseline values vector. It has the same size of the SDMcycles vector
            baselineValues = np.zeros((SDMcycles.shape[0], 1), dtype=np.float64)
            # The intercept values are not used in the code but could be used later.
            interceptDifference = np.zeros((SDMcycles.shape[0], 1), dtype=np.float64)

            outStrStuff = ''

            nfluor = rawFluor.copy()
            stupVal1 = np.zeros(spFl[0], dtype=np.float64)  # Fixme: print only
            stupVal2 = np.zeros(spFl[0], dtype=np.float64)  # Fixme: print only
            # The for loop go through all the sample matrix and make calculations well by well
            for z in range(0, spFl[0]):
                if verbose:
                    print('React: ' + str(z))
                outFFile += "----Sample: " + str(z) + "\n"  #"{0:0.6f}".format(xxelex) + "\t"
                # If there is a "no amplification" error, there is no baseline value calculated and it is automatically the
                # minimum fluorescence value assigned as baseline value for the considered reaction :
                if not vecNoAmplification[z]:
                    #  Make sure baseline is overestimated, without using slope criterion
                    #  increase baseline per cycle till eff > 2 or remaining loglin points < PointsInWoL
                    #  fastest when bgarr is directly set to 5 point below fstop
                    start = fstop[z]

                    # Find the first value that is not NaN
                    firstNotNaN = 1  # Cycles so +1 to array
                    while np.isnan(nfluor[z, firstNotNaN - 1]) and firstNotNaN < fstop[z]:
                        firstNotNaN += 1

                    subtrCount = 5
                    while subtrCount > 0 and start > firstNotNaN:
                        start -= 1
                        if not np.isnan(rawFluor[z, start - 1]):
                            subtrCount -= 1

                    bgarr[z] = 0.99 * rawFluor[z, start - 1]
                    nfluor[z] = rawFluor[z] - bgarr[z]
                    nfluor[np.isnan(nfluor)] = 0
                    nfluor[nfluor <= 0.00000001] = np.nan
                    #  baseline is now certainly too high

                    #  1. extend line downwards from fstop[] till slopelow < slopehigh of bgarr[] < sampmin[]
                    ntrials = 0
                    while True:
                        ntrials += 1
                        fstop[z] = _findLRstop(nfluor, z)
                        [fstart[z], fstart2[z]] = _findLRstart(nfluor, z, fstop[z])
                        if fstop[z] - fstart2[z] > 2:
                            [slopelow, slopehigh] = _ttestslopes(nfluor, z, fstop, fstart2)
                            defbgarr[z] = bgarr[z]
                        else:
                            break

                        if slopelow >= slopehigh:
                            bgarr[z] *= 0.99
                            nfluor[z] = rawFluor[z] - bgarr[z]
                            nfluor[np.isnan(nfluor)] = 0
                            nfluor[nfluor <= 0.00000001] = np.nan

                        if (slopelow < slopehigh or
                            bgarr[z] < 0.95 * sampmin[z] or
                            ntrials > 1000):
                            break

                    if bgarr[z] < 0.95 * sampmin[z]:
                        vecBaselineError[z] = True

                    # 2. fine tune slope of total line
                    stepval = 0.005 * bgarr[z]
                    basestep = 1.0
                    ntrials = 0
                    trialstoshift = 0
                    thisslopediff = 10
                    thissigndiff = 0
                    SlopeHasShifted = False
                    while True:
                        ntrials += 1
                        trialstoshift += 1
                        if trialstoshift > 10 and not SlopeHasShifted:
                            basestep *= 2
                            trialstoshift = 0

                        lastsigndiff = thissigndiff
                        lastslopediff = thisslopediff
                        defbgarr[z] = bgarr[z]
                        # apply baseline
                        nfluor[z] = rawFluor[z] - bgarr[z]
                        nfluor[np.isnan(nfluor)] = 0
                        nfluor[nfluor <= 0.00000001] = np.nan
                        fstop[z] = _findLRstop(nfluor, z)
                        [fstart[z], fstart2[z]] = _findLRstart(nfluor, z, fstop[z])
                        outStrStuff += "loop " + str(z + 1) + " " + str(fstop[z]) + "-" + str(fstart[z]) + " " + str(
                            nfluor[z, fstop[z] - 1]) + " " + str(nfluor[z, fstart[z] - 1]) + "\n"

                        if fstop[z] - fstart2[z] > 2:
                            [slopelow, slopehigh] = _ttestslopes(nfluor, z, fstop, fstart2)
                            thisslopediff = np.abs(slopelow - slopehigh)
                            if (slopelow - slopehigh) > 0.0:
                                thissigndiff = 1
                            else:
                                thissigndiff = -1
                            # start with baseline that is too low: lowerslope is low
                            if slopelow < slopehigh:
                                # increase baseline
                                bgarr[z] += basestep * stepval
                            else:
                                # crossed right baseline
                                # go two steps back
                                bgarr[z] -= basestep * stepval * 2
                                # decrease stepsize
                                basestep /= 2
                                SlopeHasShifted = True

                        else:
                            break

                        if (((np.abs(thisslopediff - lastslopediff) < 0.00001) and
                             (thissigndiff == lastsigndiff) and SlopeHasShifted) or
                            (np.abs(thisslopediff) < 0.0001) or
                            (ntrials > 1000)):
                            break

                    # reinstate samples that reach the slopediff criterion within 0.9 * sampmin
                    if thisslopediff < 0.0001 and defbgarr[z] > 0.9 * sampmin[z]:
                        vecBaselineError[z] = False

                    # 3: skip sample when fluor[fstop]/fluor[fstart] < 20
                    # Fixme: if RelaxLogLinLengthRG.ItemIndex = 0 then
                    if True:
                        loglinlen = 20.0
                    else:
                        loglinlen = 10.0

                    outStrStuff += "log " + str(z + 1) + " " + str(fstop[z]) + "-" + str(fstart[z]) + " " + str(nfluor[z, fstop[z] - 1]) + " " + str(nfluor[z, fstart[z] - 1]) + "\n"
                    if nfluor[z, fstop[z] - 1] / nfluor[z, fstart2[z] - 1] < loglinlen:
                        vecShortloglin[z] = True

                    pcreff[z] = np.power(10, slopehigh)

                    stupVal1[z] = nfluor[z, fstop[z] - 1]
                    stupVal2[z] = nfluor[z, fstart2[z] - 1]
                else:
                    vecSkipSample[z] = True
                    defbgarr[z] = 0.99 * sampmin[z]
                    nfluor[z] = rawFluor[z] - defbgarr[z]
                    # Fixme is this intended for all??
                    nfluor[np.isnan(nfluor)] = 0
                    nfluor[nfluor <= 0.00000001] = np.nan

                    fstop[z] = -1
                    fstart[z] = -1

                    pcreff[z] = np.nan

                    stupVal1[z] = -1.7234
                    stupVal2[z] = -1.7234

                if vecBaselineError[z]:
                    vecSkipSample[z] = True

            bgarr = defbgarr

            baselineCorrectedData = nfluor

            with open("linPas/_np_commandline.txt", "w") as f:
                f.write(outStrStuff)

            _numpyTwoAxisSave(stupVal1, "linPas/numpy_fstop_3.tsv")
            _numpyTwoAxisSave(stupVal2, "linPas/numpy_fstart_3.tsv")
            _numpyTwoAxisSave(nfluor, "linPas/numpy_fluor_data_3.tsv")
            _numpyTwoAxisSave(pcreff, "linPas/numpy_pcreff_3.tsv")
            _numpyTwoAxisSave(vecNoAmplification, "linPas/numpy_vecNoAmplification_3.tsv")
            _numpyTwoAxisSave(vecNoPlateau, "linPas/numpy_vecNoPlateau_3.tsv")
            _numpyTwoAxisSave(vecBaselineError, "linPas/numpy_vecBaselineError_3.tsv")
            _numpyTwoAxisSave(vecSkipSample, "linPas/numpy_vecSkipSample_3.tsv")
            _numpyTwoAxisSave(vecShortloglin, "linPas/numpy_vecShortloglin_3.tsv")


        if saveBaslineCorr:
            rawArr = [[header[0][rar_id], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl]]]
            for i in range(0, spFl[1]):
                rawArr[0].append(i + 1)
            for i in range(0, spFl[0]):
                rawArr.append([res[i][rar_id], res[i][rar_sample], res[i][rar_tar], res[i][rar_excl]])
                for k in range(0, spFl[1]):
                    rawArr[i + 1].append(float(baselineCorrectedData[i, k]))
            finalData["baselineCorrectedData"] = rawArr

        getmeaneff, effvar, vecIsUsedInWoL = _GetMeanEff(None, None, pcreff, vecSkipSample, vecNoPlateau, vecShortloglin, vecIsUsedInWoL)
        grpftval = upwin[0] - np.log10(getmeaneff)

        # See if cq values are stable with a modified baseline
        checkfluor = np.zeros(spFl, dtype=np.float64)
        checkstop = FDMcycles.copy()  # Fixme: only for correct length
        checkstart = FDMcycles.copy()  # Fixme: only for correct length
        checkstart2 = FDMcycles.copy()  # Fixme: only for correct length
        ctvals = np.zeros(spFl[0], dtype=np.float64)
        pcreff = np.zeros(spFl[0], dtype=np.float64)
        nnulls = np.zeros(spFl[0], dtype=np.float64)
        ninclu = np.zeros(spFl[0], dtype=np.int)
        correl = np.zeros(spFl[0], dtype=np.float64)

        for z in range(0, spFl[0]):
            if vecShortloglin[z] and not vecNoAmplification[z]:
             #    print("loglin : " + str(z))
                # No recalculation needed:
                # checkfluor[z] = rawFluor[z] - bgarr[z]
                # checkfluor[np.isnan(checkfluor)] = 0
                # checkfluor[checkfluor <= 0.00000001] = np.nan
                # checkstop[z] = _findLRstop(checkfluor, z)
                # [checkstart[z], checkstart2[z]] = _findLRstart(checkfluor, z, checkstop[z])


            #    for i in range(0, spFl[1]):
            #        if np.abs(nfluor[z, i] - checkfluor[z, i]) > 0.0000000000001:
            #            print("Size Diff" + str(z) + " : " + str(i + 1) + " " + str(nfluor[z, i]) + " " + str(checkfluor[z, i]))
            #        if not (checkstop[z] == fstop[z]):
            #            print("Stop Diff" + str(z) + " : " + str(i + 1) + " " + str(fstop[z]) + " " + str(checkstop[z]))
            #        if not (checkstop[z] == fstop[z]):
            #            print("Start Diff" + str(z) + " : " + str(i + 1) + " " + str(fstart2[z]) + " " + str(checkstart2[z]))

                uplim = np.power(10, upwin[0])  # Fixme: No logs
                lowlim = np.power(10, lowwin[0])

                checkfluor[z] = rawFluor[z] - 1.05 * bgarr[z]
                checkfluor[np.isnan(checkfluor)] = 0
                checkfluor[checkfluor <= 0.00000001] = np.nan

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    maxFlour = np.nanmax(checkfluor)

                if np.isnan(maxFlour):
                    CtShiftUp, xpcreff, xnnulls, xninclu, xcorrel = _do_cascade(nfluor, z, uplim, lowlim, grpftval)
                else:
                    CtShiftUp, xpcreff, xnulls, xninclu, xcorrel = _do_cascade(checkfluor, z, uplim, lowlim, grpftval)

                checkfluor[z] = rawFluor[z] - 0.95 * bgarr[z]
                checkfluor[np.isnan(checkfluor)] = 0
                checkfluor[checkfluor <= 0.00000001] = np.nan
                CtShiftDown, xpcreff, xnnulls, xninclu, xcorrel = _do_cascade(checkfluor, z, uplim, lowlim, grpftval)

                if np.abs(CtShiftUp - CtShiftDown) > 1.0:
                    vecBaselineError[j] = True
                    vecSkipSample[j] = True
                else:
                    if not vecBaselineError[j]:
                        vecSkipSample[j] = False

        vecSkipSample[vecExcludedByUser] = True
        # Update the window
        meanmax = _GetMeanMax(nfluor, None, vecSkipSample, vecNoPlateau)
        upwin[0] = np.log10(0.1 * meanmax)
        lowwin[0] = np.log10(0.1 * meanmax / 16.0)
        logfluor = np.log10(nfluor)
        yaxismax = np.nanmax(logfluor)
        yaxismin = np.nanmin(logfluor)
        if yaxismin < yaxismax - 5:
            yaxismin = yaxismax - 5

        # CheckNoisiness
        skipgroup = False
        maxlim = _GetMeanFluStop(nfluor, None, None, fstop, vecSkipSample, vecNoPlateau)
        if maxlim > 0.0:
            maxlim = np.log10(maxlim)
        else:
            skipgroup = True
        thismeaneff = 1.0

        PointsInWoL = 4

        if not skipgroup:
            step = PointsInWoL * _GetLogStepStop(nfluor, None, None, fstop, vecSkipSample, vecNoPlateau)
            upwin, lowwin = _ApplyLogWindow(None, maxlim, step, upwin, lowwin, yaxismax, yaxismin)
            # Fixme: no trunc needed
            grpftval = np.log10(0.5 * np.trunc(1000 * np.power(10, upwin[0])) / 1000)
            tctvals, tpcreff, tnnulls, tninclu, tcorrel = _do_all_cascade(nfluor, None, None, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin,
                                                                          grpftval, vecSkipSample, vecNoPlateau)
            thismeaneff, effvar, vecIsUsedInWoL = _GetMeanEff(None, None, tpcreff, vecSkipSample, vecNoPlateau, vecShortloglin, vecIsUsedInWoL)
            if thismeaneff < 1.001:
                skipgroup = True

        SetIsNoisy = False
        if not skipgroup:
            foldwidth = np.log10(np.power(thismeaneff, PointsInWoL))
            upwin, lowwin = _ApplyLogWindow(None, maxlim, foldwidth, upwin, lowwin, yaxismax, yaxismin)

        # FixMe: Still a lot missing for noisy

        # Calculate the WoL per group
        for t in range(1, targetsCount):
            upwin[t] = upwin[0]
            lowwin[t] = lowwin[0]

        for t in range(1, targetsCount):
            ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin, vecIsUsedInWoL = _Set_WoL(nfluor, t, vecTarget, PointsInWoL, ctvals, pcreff, nnulls, ninclu, correl, upwin, lowwin, yaxismax, yaxismin, fstop, fstart2, grpftval, vecSkipSample, vecNoPlateau, vecShortloglin, vecIsUsedInWoL)
            # AssignNoPlateau(k)

        # Assign no plateau




        print("----------------------------")


        _numpyTwoAxisSave(vecIsUsedInWoL, "linPas/numpy_vecIsUsedInWoL_3.tsv")


        with open("matlab/numpy_variables.txt", "w") as f:
            f.write(outFFile)

        # Fixme: Delete
        stop_time = dt.datetime.now() - start_time
        print("Done Baseline: " + str(stop_time) + "sec")

        #################################
        # CALCULATION OF A WOL PER GENE #
        #################################
        # Fixme: update cq in RDML

        for i in range(0, len(res)):
           # if not vecExcludedByUser[i]:
            res[i][rar_baseline] = np.log10(bgarr[i])
            res[i][rar_lower_limit] = np.power(10, lowwin[vecTarget[i]])
            res[i][rar_threshold] = grpftval
            res[i][rar_upper_limit] = np.power(10, upwin[vecTarget[i]])
            res[i][rar_n_included] = ninclu[i]
            res[i][rar_Cq] = ctvals[i]
            res[i][rar_indiv_PCR_eff] = pcreff[i]
            res[i][rar_N0_indiv_eff] = nnulls[i]
            res[i][rar_R2] = correl[i] * correl[i]
            res[i][rar_mean_PCR_eff] = 0  # calculThreshold
            res[i][rar_N0_mean_eff] = 0  # calculThreshold

            res[i][rar_no_amplification] = vecNoAmplification[i]
            res[i][rar_baseline_error] = vecBaselineError[i]
    #        res[i][rar_no_plateau] = vecExcludedNoPlateau[i]
    #        res[i][rar_noisy_sample] = vecExcludedNoisySample[i]
    #        res[i][rar_PCR_efficiency_outside] = vecExcludedPCREfficiency[i]

        # Fixme: Delete
        stop_time = dt.datetime.now() - start_time
        print("Done All: " + str(stop_time) + "sec")


        # optimalLowerLimit, optimalUpperLimit,



     #   print(str(calculThreshold) + " " + str(meanEfficiency) + " " + str(CVmin) + " " + str(optimalUpperLimit))
     #   print(str(optimalLowerLimit) + " " + str(optimalLowerLimit) + " " + str(optimalLowerLimit) + " " + str(optimalLowerLimit))

   #     _numpyTwoAxisSave(indexOutliers, "res_arr_5.tsv")
  #      _numpyTwoAxisSave(vecExcludedPCREfficiency, "res_arr_6.tsv")
  #      _numpyTwoAxisSave(vecExcludedNoisySample, "res_arr_7.tsv")


       # allGenes   rowRes[][rar_tar]
      #  vecExcludedByUser = [rowRes[rar_excl] for rowRes in res]
     #   rar_excl

        # print(geneNames)

        if saveResultsCSV:
            retCSV = ""
            for j in range(0, len(header[0])):
                retCSV += header[0][j] + "\t"
            retCSV = re.sub(r"\t$", "\n", retCSV)

            for i in range(0, len(res)):
                for j in range(0, len(res[i])):
                    if j in [rar_no_amplification, rar_baseline_error, rar_no_plateau, rar_noisy_sample, rar_PCR_efficiency_outside]:
                        if res[i][j]:
                            retCSV += str(res[i][j]) + "\t"
                        else:
                            retCSV += "\t"
                    else:
                        retCSV += str(res[i][j]) + "\t"
                retCSV = re.sub(r"\t$", "\n", retCSV)
            finalData["resultsCSV"] = retCSV

        if saveResultsList:
            finalData["resultsList"] = res

        return finalData


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The command line interface to the RDML-Python library.')
    parser.add_argument('-v', '--validate', metavar="data.rdml", dest='validate', help='validate file against schema')
    parser.add_argument("-d", "--doooo", dest="doooo", help="just do stuff")

    args = parser.parse_args()

    # Validate RDML file
    if args.validate:
        instValidate = Rdml()
        resValidate = instValidate.validate(filename=args.validate)
        print(resValidate)
        sys.exit(0)

    # Tryout things
    if args.doooo:
        print('Test LinRegPCR')
        rt = Rdml(args.doooo)
        xxexp = rt.experiments()
        xxrun = xxexp[0].runs()
        xxres = xxrun[0].linRegPCR(baselineCorr=True, saveRaw=True, saveBaslineCorr=True, saveResultsCSV=True, perTarget=False, verbose=False)
        if "rawData" in xxres:
            with open("test/temp_out_raw_data.tsv", "w") as f:
                xxxResStr = ""
                for xxxrow in xxres["rawData"]:
                    for xxelex in xxxrow:
                        if type(xxelex) is float:
                            xxxResStr += "{0:0.3f}".format(xxelex) + "\t"
                        else:
                            xxxResStr += str(xxelex) + "\t"
                    xxxResStr = re.sub(r"\t$", "\n", xxxResStr)
                f.write(xxxResStr)
        if "baselineCorrectedData" in xxres:
            with open("test/temp_out_baseline_corrected_data.tsv", "w") as f:
                xxxResStr = ""
                for xxxrow in xxres["baselineCorrectedData"]:
                    for xxelex in xxxrow:
                        if type(xxelex) is float:
                            xxxResStr += "{0:0.12f}".format(xxelex) + "\t"
                        else:
                            xxxResStr += str(xxelex) + "\t"
                    xxxResStr = re.sub(r"\t$", "\n", xxxResStr)
                f.write(xxxResStr)
        if "resultsCSV" in xxres:
            with open("test/temp_out_results.tsv", "w") as f:
                f.write(xxres["resultsCSV"])

        # rt.save('new.rdml')
        sys.exit(0)