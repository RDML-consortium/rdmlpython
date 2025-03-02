#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import sys
import io
import os
import re
import datetime
import zipfile
import tempfile
import argparse
import math
import warnings
import json
import csv
import numpy as np
import scipy.stats as scp
from matplotlib.backends.backend_svg import FigureCanvasSVG as figSVG
from matplotlib.figure import Figure as plt_fig
from lxml import etree as et


def get_rdml_lib_version():
    """Return the version string of the RDML library.

    Returns:
        The version string of the RDML library.
    """

    return "3.0.0"


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)


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


def _add_err(err, sep, text):
    """Get the var as float or "".

    Args:
        err: The error string to append to.
        sep: The seperator used if the string is not empty.
        text: The text to append.

    Returns:
        The error string.
    """

    if err == "":
        return text
    else:
        return err + sep + text


def _save_float(var):
    """Get the var as float or "".

    Args:
        var: The value as float.

    Returns:
        The float or "" on error.
    """

    if var == "":
        return ""
    try:
        vFloat = float(var)
    except ValueError:
        return ""
    else:
        if math.isfinite(vFloat):
            return vFloat
        else:
            return ""


def _get_copies_threshold():
    """Get the number of copies at threshold (rule of thumb).

    Returns:
        The number of copies at threshold as float.
    """

    copies = 10 * pow(1.9, 35.0) * 20.0 / 50.0
    return copies


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

    return et.Element(tag, id=id)


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
        et.SubElement(base, tag).text = text
    else:
        if text is not None and text != "":
            et.SubElement(base, tag).text = text


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
            new_node = et.Element(tag)
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
        new_node = et.Element(tag)
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
        try:
            with zipfile.ZipFile(rdmlName, 'r') as RDMLin:
                for item in RDMLin.infolist():
                    if item.filename == fileName:
                        needRewrite = True
        except zipfile.BadZipFile as e:
            needRewrite = False

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


def _niceQuantityType(txt):
    if txt == "cop":
        return "copies per microliter"
    if txt == "fold":
        return "fold change"
    if txt == "dil":
        return "dilution"
    if txt == "nMol":
        return "nanomol per microliter"
    if txt == "ng":
        return "nanogram per microliter"
    if txt == "other":
        return "other unit"
    return txt


def _niceQuantityTypeReact(txt):
    if txt == "cop":
        return "copies per reaction"
    if txt == "fold":
        return "fold change"
    if txt == "dil":
        return "dilution"
    if txt == "nMol":
        return "nanomol per reaction"
    if txt == "ng":
        return "nanogram per reaction"
    if txt == "other":
        return "other unit"
    return txt


def _wellRangeToList(wells, columns):
    """Translates a range like "B2-C3" to list ["B2","B3","C2","C3"].

    Args:
        wells: The range like "B2-C3"
        columns: Number of colums in a plate

    Returns:
        returns a list like ["B2","B3","C2","C3"].
    """

    res = []
    wellList = []
    re_found_comp = re.search(r"(\D)(\d+)-(\D)(\d+)\s*", wells)
    if re_found_comp:
        letter_s = ord(re_found_comp.group(1).upper())
        letter_e = ord(re_found_comp.group(3).upper())
        numb_s = int(re_found_comp.group(2))
        numb_e = int(re_found_comp.group(4))
        if letter_e < letter_s:
            temp_num = letter_s
            letter_s = letter_e
            letter_e = temp_num
        if numb_e < numb_s:
            temp_num = numb_s
            numb_s = numb_e
            numb_e = temp_num
        for letter in range(letter_s, letter_e + 1):
            for numb in range(numb_s, numb_e + 1):
                wellList.append(chr(letter) + str(numb))
    else:
        re_found_simp = re.search(r"(\D)(\d+)\s*", wells)
        if re_found_simp:
            wellList.append(re_found_simp.group(1) + re_found_simp.group(2))

    for well in wellList:
        re_found_simp = re.search(r"(\D)(\d+)\s*", well)
        if re_found_simp:
            letter = ord(re_found_simp.group(1).upper()) - ord("A")
            numb = int(re_found_simp.group(2))
            fin_id = numb + letter * int(columns)
            res.append(str(fin_id))
    return res


def _sampleTypeToDics(rootEl):
    """Translates all ele in samples to dic {sample Id}{target Id} with ele value.

    Args:
        rootEl: The root node of the rdml file

    Returns:
        returns a dictionary {sample Id}{target Id} with ele value.
    """
    tar = []
    ret = {}

    # Get all targets
    targets = _get_all_children(rootEl, "target")
    for node in targets:
        if 'id' in node.attrib:
            tar.append(node.attrib['id'])

    samples = _get_all_children(rootEl, "sample")
    for node in samples:
        if 'id' in node.attrib:
            currSam = node.attrib['id']
            ret[currSam] = {}
            for currTar in tar:
                ret[currSam][currTar] = "unkn"
            subEles = _get_all_children(node, "type")
            for subNode in subEles:
                if 'targetId' not in subNode.attrib:
                    for currTar in tar:
                        ret[currSam][currTar] = subNode.text
            for subNode in subEles:
                if 'targetId' in subNode.attrib:
                    ret[currSam][subNode.attrib['targetId']] = subNode.text
    return ret


def _lrp_linReg(adp_cyc_min, fluor):
    """A function which calculates the slope or the intercept by linear regression.

    Args:
        xIn: The numpy array of the cycles
        yUse: The numpy array that contains the fluorescence

    Returns:
        An array with the slope and intercept.
    """

    # Create a matrix with the cycle for each rawFluor value
    spFl = fluor.shape
    cycles = np.arange(adp_cyc_min, (adp_cyc_min + spFl[1]))
    cyclesSqared = cycles * cycles
    cyclesSum = np.sum(cycles)
    cyclesSqaredSum = np.sum(cyclesSqared)

    cycFluor = cycles * fluor
    sumFluor = np.nansum(fluor, axis=1)
    sumCycFluor = np.nansum(cycFluor, axis=1)

    ssx = cyclesSqaredSum - (cyclesSum * cyclesSum) / spFl[1]
    sxy = sumCycFluor - (cyclesSum * sumFluor) / spFl[1]

    return sxy / ssx
    # slope = sxy / ssx
    # intercept = (sumFluor / n) - slope * (sumCyc / n)
    # return [slope, intercept]


def _lrp_findStopCyc(fluor, aRow):
    """Find the stop cycle of the log lin phase in fluor.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on

    Returns:
        An int with the stop cycle.
    """

    # Take care of nan values
    validTwoLessCyc = 3  # Cycles so +1 to array
    while (validTwoLessCyc <= fluor.shape[1] and
           (np.isnan(fluor[aRow, validTwoLessCyc - 1]) or
            np.isnan(fluor[aRow, validTwoLessCyc - 2]) or
            np.isnan(fluor[aRow, validTwoLessCyc - 3]))):
        validTwoLessCyc += 1

    # First and Second Derivative values calculation
    fluorShift = np.roll(fluor[aRow], 1, axis=0)  # Shift to right - real position is -0.5
    fluorShift[0] = np.nan
    firstDerivative = fluor[aRow] - fluorShift
    if np.isfinite(firstDerivative).any():
        FDMaxCyc = np.nanargmax(firstDerivative, axis=0) + 1  # Cycles so +1 to array
    else:
        return fluor.shape[1]

    firstDerivativeShift = np.roll(firstDerivative, -1, axis=0)  # Shift to left
    firstDerivativeShift[-1] = np.nan
    secondDerivative = firstDerivativeShift - firstDerivative

    if FDMaxCyc + 2 <= fluor.shape[1]:
        # Only add two cycles if there is an increase without nan
        if (not np.isnan(fluor[aRow, FDMaxCyc - 1]) and
                not np.isnan(fluor[aRow, FDMaxCyc]) and
                not np.isnan(fluor[aRow, FDMaxCyc + 1]) and
                fluor[aRow, FDMaxCyc + 1] > fluor[aRow, FDMaxCyc] > fluor[aRow, FDMaxCyc - 1]):
            FDMaxCyc += 2
    else:
        FDMaxCyc = fluor.shape[1]

    maxMeanSD = 0.0
    stopCyc = fluor.shape[1]

    for cycInRange in range(validTwoLessCyc, FDMaxCyc):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            tempMeanSD = np.mean(secondDerivative[cycInRange - 2: cycInRange + 1], axis=0)
        # The > 0.000000000001 is to avoid float differences to the pascal version
        if not np.isnan(tempMeanSD) and (tempMeanSD - maxMeanSD) > 0.000000000001:
            maxMeanSD = tempMeanSD
            stopCyc = cycInRange
    if stopCyc + 2 >= fluor.shape[1]:
        stopCyc = fluor.shape[1]

    return stopCyc


def _lrp_findStartCyc(fluor, row, stopCyc):
    """A function which finds the start cycle of the log lin phase in fluor.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on
        stopCyc: The stop cycle

    Returns:
        An array [int, int] with the start cycle and the fixed start cycle.
    """

    startCyc = stopCyc - 1
    if startCyc < 1 or np.isnan(fluor[row, startCyc - 1]):
        return [startCyc, startCyc]

    # As long as there are no NaN and new values are increasing
    while (startCyc > 1 and
           stopCyc - startCyc < 11 and
           not np.isnan(fluor[row, startCyc - 2]) and
           fluor[row, startCyc - 2] < fluor[row, startCyc - 1]):
        startCyc -= 1

    startCycFix = startCyc
    startStep = np.log10(fluor[row, startCyc]) - np.log10(fluor[row, startCyc - 1])
    stopStep = np.log10(fluor[row, stopCyc - 1]) - np.log10(fluor[row, stopCyc - 2])
    if startStep > 1.1 * stopStep:
        startCycFix += 1

    return [startCyc, startCycFix]


def _lrp_findStartCyc2(fluor, row, stopCyc):
    """A function which finds the start cycle of the log lin phase in fluor.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on
        stopCyc: The stop cycle

    Returns:
        An array [int, int] with the start cycle and the fixed start cycle.
    """

    startCyc = stopCyc
    if (startCyc < 1 or 
        np.isnan(fluor[row, startCyc - 1]) or
        fluor[row, startCyc - 1] <= 0.0):
        return [startCyc, startCyc]

    # As long as there are no NaN and new values are increasing
    while (startCyc > 1 and
           stopCyc - startCyc < 11 and  # TODO set to 6!
           not np.isnan(fluor[row, startCyc - 2]) and
           fluor[row, startCyc - 2] > 0.0 and
           fluor[row, startCyc - 2] < fluor[row, startCyc - 1]):
        startCyc -= 1

    startCycFix = startCyc
    if startCyc < stopCyc:
        startStep = np.log10(fluor[row, startCyc]) - np.log10(fluor[row, startCyc - 1])
        stopStep = np.log10(fluor[row, stopCyc - 1]) - np.log10(fluor[row, stopCyc - 2])
        if startStep > 1.1 * stopStep:
            startCycFix += 1

    return [startCyc, startCycFix]


def _lrp_testSlopes(fluor, aRow, stopCyc, startCycFix):
    """Splits the values and calculates a slope for the upper and the lower half.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on
        stopCyc: The stop cycle
        startCycFix: The start cycle

    Returns:
        An array with [slopelow, slopehigh].
    """

    # Both start with full range
    center = (stopCyc[aRow] + startCycFix[aRow]) / 2.0
    loopStart = [startCycFix[aRow], int(center + 0.6)]
    loopStop = [int(center), stopCyc[aRow]]

    # basic regression per group
    ssx = [0, 0]
    sxy = [0, 0]
    slope = [0, 0]
    for j in range(0, 2):
        sumx = 0.0
        sumy = 0.0
        sumx2 = 0.0
        sumxy = 0.0
        nincl = 0.0
        for i in range(loopStart[j], loopStop[j] + 1):
            sumx += i
            sumy += np.log10(fluor[aRow, i - 1])
            sumx2 += i * i
            sumxy += i * np.log10(fluor[aRow, i - 1])
            nincl += 1
        if nincl != 0.0:
            ssx[j] = sumx2 - sumx * sumx / nincl
            sxy[j] = sumxy - sumx * sumy / nincl
            slope[j] = sxy[j] / ssx[j]
        else:
            slope[j] = -999.9

    return [slope[0], slope[1]]


def _lrp_testSlopes2(fluor):
    """Splits the values and calculates a slope for the upper and the lower half.

    Args:
        fluor: The array with the fluorescence values

    Returns:
        An array with [slopelow, slopehigh].
    """

    # Split the list
    center = len(fluor) / 2.0
    loopStart = [0, int(center)]
    loopStop = [int(center - 0.4), len(fluor) - 1]

    # basic regression per section
    ssx = [0, 0]
    sxy = [0, 0]
    slope = [0, 0]
    for j in range(0, 2):
        sumx = 0.0
        sumy = 0.0
        sumx2 = 0.0
        sumxy = 0.0
        nincl = 0.0
        for i in range(loopStart[j], loopStop[j] + 1):
            sumx += i
            sumy += np.log10(fluor[i])
            sumx2 += i * i
            sumxy += i * np.log10(fluor[i])
            nincl += 1
        if nincl != 0.0:
            ssx[j] = sumx2 - sumx * sumx / nincl
            sxy[j] = sumxy - sumx * sumy / nincl
            slope[j] = sxy[j] / ssx[j]
        else:
            slope[j] = -999.9

    return [slope[0], slope[1]]


def _lrp_testSlopes3(fluor):
    """Splits the values and calculates a slope for the upper and the lower half.

    Args:
        fluor: The array with the fluorescence values

    Returns:
        An array with [slopelow, slopehigh].
    """

    # Split the list
    center = len(fluor) / 2.0
    loopStart = 0
    loopStop = len(fluor)

    # basic regression per section
    ssx = 0.0
    sxy = 0.0
    slope = 0.0

    sumx = 0.0
    sumy = 0.0
    sumx2 = 0.0
    sumxy = 0.0
    nincl = 0.0
    for i in range(loopStart, loopStop):
        sumx += i
        sumy += np.log10(fluor[i])
        sumx2 += i * i
        sumxy += i * np.log10(fluor[i])
        nincl += 1
    if nincl != 0.0:
        ssx = sumx2 - sumx * sumx / nincl
        sxy = sumxy - sumx * sumy / nincl
        slope = sxy / ssx
    else:
        slope = -999.9

    return slope





def _lrp_lastCycMeanMax(fluor, vecSkipSample, vecNoPlateau):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        fluor: The array with the fluorescence values
        vecSkipSample: Skip the sample
        vecNoPlateau: Sample has no plateau

    Returns:
        An float with the max mean.
    """

    # Ignore all nan slices, to fix them below
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        maxFlour = np.nanmax(fluor[:, -11:], axis=1)
        maxFlour[vecSkipSample] = np.nan
        maxFlour[vecNoPlateau] = np.nan
        maxMean = np.nanmean(maxFlour)
    if np.isnan(maxMean):
        maxMean = np.nanmax(maxFlour)

    return maxMean


def _lrp_meanPcrEff(tarGroup, vecTarget, pcrEff, vecSkipSample, vecNoPlateau, vecShortLogLin):
    """A function which calculates the mean efficiency of the selected target group excluding bad ones.

    Args:
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        pcrEff: The array with the PCR efficiencies
        vecSkipSample: Skip the sample
        vecNoPlateau: True if there is no plateau
        vecShortLogLin: True indicates a short log lin phase

    Returns:
        An array with [meanPcrEff, pcrEffVar].
    """

    cnt = 0
    sumEff = 0.0
    sumEff2 = 0.0
    for j in range(0, len(pcrEff)):
        if tarGroup is None or tarGroup == vecTarget[j]:
            if (not (vecSkipSample[j] or vecNoPlateau[j] or vecShortLogLin[j])) and pcrEff[j] > 1.0:
                cnt += 1
                sumEff += pcrEff[j]
                sumEff2 += pcrEff[j] * pcrEff[j]

    if cnt > 1:
        meanPcrEff = sumEff / cnt
        pcrEffVar = (sumEff2 - (sumEff * sumEff) / cnt) / (cnt - 1)
    else:
        meanPcrEff = 1.0
        pcrEffVar = 100

    return [meanPcrEff, pcrEffVar]


def _lrp_startStopInWindow(fluor, row, startCyc, startCycFix, stopCyc, upWin, lowWin):
    """Find the start and the stop of the part of the curve which is inside the window.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on
        startCyc: The start cycle of log phase
        startCycFix: The fixed start cycle of log phase
        upWin: The upper limit of the window
        lowWin: The lower limit of the window

    Returns:
        The int startWinCyc, stopWinCyc and the bool notInWindow.
    """

    startWinCyc = 0
    stopWinCyc = 0

    if np.isfinite(fluor[row, startCycFix[row] - 1:]).any():
        # stopMaxCyc = np.nanargmax(fluor[row, startCycFix[row] - 1:]) + startCycFix[row]
        # TODO: Why expanding to max?
        stopMaxCyc = stopCyc[row]
        #if stopCyc[row] != stopMaxCyc:
          #  print("STOP issue: " + str(stopCyc[row]) + " vs " + str(stopMaxCyc))
    else:
        return startCyc[row], startCyc[row], True

    # If is true if outside the window
    if fluor[row, startCyc[row] - 1] > upWin or fluor[row, stopMaxCyc - 1] < lowWin:
        notInWindow = True
        if fluor[row, startCyc[row] - 1] > upWin:
            startWinCyc = startCyc[row]
            stopWinCyc = startCyc[row]
        if fluor[row, stopMaxCyc - 1] < lowWin:
            startWinCyc = stopMaxCyc
            stopWinCyc = stopMaxCyc
    else:
        notInWindow = False
        # look for stopWinCyc
        if fluor[row, stopMaxCyc - 1] < upWin:
            stopWinCyc = stopMaxCyc
        else:
            for i in range(stopMaxCyc, startCyc[row], -1):
                if fluor[row, i - 1] > upWin > fluor[row, i - 2]:
                    stopWinCyc = i - 1
        # look for startWinCyc
        if fluor[row, startCycFix[row] - 1] > lowWin:
            startWinCyc = startCycFix[row]
        else:
            for i in range(stopMaxCyc, startCyc[row], -1):
                if fluor[row, i - 1] > lowWin > fluor[row, i - 2]:
                    startWinCyc = i
    return startWinCyc, stopWinCyc, notInWindow


def _lrp_paramInWindow(fluor, aRow, startCyc, startCycFix, stopCyc, upWin, lowWin):
    """Calculates slope, nNull, PCR efficiency and mean x/y for the curve part in the window.

    Args:
        fluor: The array with the fluorescence values
        aRow: The row to work on
        startCyc: The start cycle of log phase
        startCycFix: The fixed start cycle of log phase
        stopCyc: The stop cycle of log phase
        upWin: The upper limit of the window
        lowWin: The lower limit of the window

    Returns:
        The calculated values: indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl.
    """

    startWinCyc, stopWinCyc, notInWindow = _lrp_startStopInWindow(fluor, aRow, startCyc, startCycFix, stopCyc, upWin, lowWin)

    sumx = 0.0
    sumy = 0.0
    sumx2 = 0.0
    sumy2 = 0.0
    sumxy = 0.0
    nincl = 0.0
    ssx = 0.0
    ssy = 0.0
    sxy = 0.0
    for i in range(startWinCyc, stopWinCyc + 1):
        fluorSamp = fluor[aRow, i - 1]
        if not np.isnan(fluorSamp):
            logFluorSamp = np.log10(fluorSamp)
            sumx += i
            sumy += logFluorSamp
            sumx2 += i * i
            sumy2 += logFluorSamp * logFluorSamp
            sumxy += i * logFluorSamp
            nincl += 1

    if nincl > 1:
        ssx = sumx2 - sumx * sumx / nincl
        ssy = sumy2 - sumy * sumy / nincl
        sxy = sumxy - sumx * sumy / nincl

    if ssx > 0.0 and ssy > 0.0 and nincl > 0.0:
        cslope = sxy / ssx
        cinterc = sumy / nincl - cslope * sumx / nincl
        correl = sxy / np.sqrt(ssx * ssy)
        indMeanX = sumx / nincl
        indMeanY = sumy / nincl
        pcrEff = np.power(10, cslope)
        nnulls = np.power(10, cinterc)
    else:
        correl = np.nan
        indMeanX = np.nan
        indMeanY = np.nan
        pcrEff = np.nan
        nnulls = np.nan

    if notInWindow:
        ninclu = 0
    else:
        ninclu = stopWinCyc - startWinCyc + 1

    return indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl


def _lrp_allParamInWindow(fluor, tarGroup, vecTarget, indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl, startCyc, startCycFix, stopCyc, upWin, lowWin, vecNoAmplification, vecBaselineError):
    """A function which calculates the mean of the max fluor in the last ten cycles.

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        indMeanX: The vector with the x mean position
        indMeanY: The vector with the y mean position
        pcrEff: The array with the PCR efficiencies
        nnulls: The array with the calculated nnulls
        ninclu: The array with the calculated ninclu
        correl: The array with the calculated correl
        startCyc: The start cycle of log phase
        startCycFix: The fixed start cycle of log phase
        stopCyc: The stop cycle of log phase
        upWin: The upper limit of the window
        lowWin: The lower limit of the window
        vecNoAmplification: True if there is a amplification error
        vecBaselineError: True if there is a baseline error

    Returns:
        An array with [indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl].
    """

    for row in range(0, fluor.shape[0]):
        if tarGroup is None or tarGroup == vecTarget[row]:
            if not (vecNoAmplification[row] or vecBaselineError[row]):
                if tarGroup is None:
                    indMeanX[row], indMeanY[row], pcrEff[row], nnulls[row], ninclu[row], correl[row] = _lrp_paramInWindow(fluor, row, startCyc, startCycFix, stopCyc, upWin[0], lowWin[0])
                else:
                    indMeanX[row], indMeanY[row], pcrEff[row], nnulls[row], ninclu[row], correl[row] = _lrp_paramInWindow(fluor, row, startCyc, startCycFix, stopCyc, upWin[tarGroup], lowWin[tarGroup])
            else:
                correl[row] = np.nan
                indMeanX[row] = np.nan
                indMeanY[row] = np.nan
                pcrEff[row] = np.nan
                nnulls[row] = np.nan
                ninclu[row] = 0

    return indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl


def _lrp_meanStopFluor(fluor, tarGroup, vecTarget, stopCyc, vecSkipSample, vecNoPlateau):
    """Return the mean of the stop fluor or the max fluor if all rows have no plateau.

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        stopCyc: The vector with the stop cycle of the log lin phase
        vecSkipSample: Skip the sample
        vecNoPlateau: True if there is no plateau

    Returns:
        The meanMax fluorescence.
    """

    meanMax = 0.0
    maxFluor = 0.0000001
    cnt = 0
    if tarGroup is None:
        for aRow in range(0, fluor.shape[0]):
            if not vecSkipSample[aRow]:
                if not vecNoPlateau[aRow]:
                    cnt += 1
                    meanMax += fluor[aRow, stopCyc[aRow] - 1]
                else:
                    for i in range(0, fluor.shape[1]):
                        if fluor[aRow, i] > maxFluor:
                            maxFluor = fluor[aRow, i]
    else:
        for aRow in range(0, fluor.shape[0]):
            if tarGroup == vecTarget[aRow] and not vecSkipSample[aRow]:
                if not vecNoPlateau[aRow]:
                    cnt += 1
                    meanMax += fluor[aRow, stopCyc[aRow] - 1]
                else:
                    for i in range(0, fluor.shape[1]):
                        if fluor[aRow, i] > maxFluor:
                            maxFluor = fluor[aRow, i]

    if cnt > 0:
        meanMax = meanMax / cnt
    else:
        meanMax = maxFluor
    return meanMax


def _lrp_maxStartFluor(fluor, tarGroup, vecTarget, startCyc, vecSkipSample):
    """Return the maximum of the start fluorescence

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        startCyc: The vector with the start cycle of the log lin phase
        vecSkipSample: Skip the sample

    Returns:
        The maxStart fluorescence.
    """
    maxStart = -10.0
    if tarGroup is None:
        for aRow in range(0, fluor.shape[0]):
            if not vecSkipSample[aRow]:
                if fluor[aRow, startCyc[aRow] - 1] > maxStart:
                    maxStart = fluor[aRow, startCyc[aRow] - 1]
    else:
        for aRow in range(0, fluor.shape[0]):
            if tarGroup == vecTarget[aRow] and not vecSkipSample[aRow]:
                if fluor[aRow, startCyc[aRow] - 1] > maxStart:
                    maxStart = fluor[aRow, startCyc[aRow] - 1]

    return 0.999 * maxStart


def _lrp_setLogWin(tarGroup, newUpWin, foldWidth, upWin, lowWin, maxFluorTotal, minFluorTotal):
    """Sets a new window and ensures its within the total fluorescence values.

    Args:
        tarGroup: The target number
        newUpWin: The new upper window
        foldWidth: The foldWith to the lower window
        upWin: The upper window fluorescence
        lowWin: The lower window fluorescence
        maxFluorTotal: The maximum fluorescence over all rows
        minFluorTotal: The minimum fluorescence over all rows

    Returns:
        The arrays upWin and lowWin.
    """
    # No rounding needed, only present for exact identical output with Pascal version
    tempUpWin = np.power(10, np.round(1000 * newUpWin) / 1000)
    tempLowWin = np.power(10, np.round(1000 * (newUpWin - foldWidth)) / 1000)

    tempUpWin = np.minimum(tempUpWin, maxFluorTotal)
    tempUpWin = np.maximum(tempUpWin, minFluorTotal)
    tempLowWin = np.minimum(tempLowWin, maxFluorTotal)
    tempLowWin = np.maximum(tempLowWin, minFluorTotal)

    if tarGroup is None:
        upWin[0] = tempUpWin
        lowWin[0] = tempLowWin
    else:
        upWin[tarGroup] = tempUpWin
        lowWin[tarGroup] = tempLowWin

    return upWin, lowWin


def _lrp_logStepStop(fluor, tarGroup, vecTarget, stopCyc, vecSkipSample, vecNoPlateau):
    """Calculates the log of the fluorescence increase at the stop cycle.

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        stopCyc: The vector with the stop cycle of the log lin phase
        vecSkipSample: True if row should be skipped
        vecNoPlateau: True if there is no plateau

    Returns:
        An array with [indMeanX, indMeanY, pcrEff, nnulls, ninclu, correl].
    """
    cnt = 0
    step = 0.0
    for aRow in range(0, fluor.shape[0]):
        if (tarGroup is None or tarGroup == vecTarget[aRow]) and not (vecSkipSample[aRow] or vecNoPlateau[aRow]):
            cnt += 1
            step += np.log10(fluor[aRow, stopCyc[aRow] - 1]) - np.log10(fluor[aRow, stopCyc[aRow] - 2])
    if cnt > 0:
        step = step / cnt
    else:
        step = np.log10(1.8)
    return step


def _lrp_setWoL(fluor, tarGroup, vecTarget, pointsInWoL, indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl,
                upWin, lowWin, maxFluorTotal, minFluorTotal, stopCyc, startCyc, startCycFix, threshold,
                vecNoAmplification, vecBaselineError, vecSkipSample, vecNoPlateau, vecShortLogLin, vecIsUsedInWoL):
    """Find the window with the lowest variation in PCR efficiency and calculate its values.

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        pointsInWoL: The number of points in the window
        indMeanX: The vector with the x mean position
        indMeanY: The vector with the y mean position
        pcrEff: The array with the PCR efficiencies
        nNulls: The array with the calculated nNulls
        nInclu: The array with the calculated nInclu
        correl: The array with the calculated correl
        upWin: The upper limit of the window
        lowWin: The lower limit of the window
        maxFluorTotal: The maximum fluorescence over all rows
        minFluorTotal: The minimum fluorescence over all rows
        stopCyc: The vector with the stop cycle of the log lin phase
        startCyc: The vector with the start cycle of the log lin phase
        threshold: The threshold fluorescence
        vecNoAmplification: True if there is a amplification error
        vecBaselineError: True if there is a baseline error
        vecSkipSample: Skip the sample
        vecNoPlateau: True if there is no plateau
        vecShortLogLin: True indicates a short log lin phase
        vecIsUsedInWoL: True if used in the WoL

    Returns:
        The values indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL.
    """
    skipGroup = False
    stepSize = 0.2  # was 0.5, smaller steps help in finding WoL
    # Keep 60 calculated results
    memVarEff = np.zeros(60, dtype=np.float64)
    memUpWin = np.zeros(60, dtype=np.float64)
    memFoldWidth = np.zeros(60, dtype=np.float64)

    maxFluorWin = _lrp_meanStopFluor(fluor, tarGroup, vecTarget, stopCyc, vecSkipSample, vecNoPlateau)
    if maxFluorWin > 0.0:
        maxFluorWin = np.log10(maxFluorWin)
    else:
        skipGroup = True
    minFluorLim = _lrp_maxStartFluor(fluor, tarGroup, vecTarget, startCyc, vecSkipSample)
    if minFluorLim > 0.0:
        minFluorLim = np.log10(minFluorLim)
    else:
        skipGroup = True

    checkMeanEff = 1.0
    if not skipGroup:
        foldWidth = pointsInWoL * _lrp_logStepStop(fluor, tarGroup, vecTarget, stopCyc, vecSkipSample, vecNoPlateau)
        upWin, lowWin = _lrp_setLogWin(tarGroup, maxFluorWin, foldWidth, upWin, lowWin, maxFluorTotal, minFluorTotal)

        _unused, _unused2, checkPcrEff, _unused3, _unused4, _unused5 = _lrp_allParamInWindow(fluor, tarGroup, vecTarget,
                                                                                             indMeanX, indMeanY, pcrEff,
                                                                                             nNulls, nInclu, correl,
                                                                                             startCyc, startCycFix, stopCyc,
                                                                                             upWin, lowWin,
                                                                                             vecNoAmplification,
                                                                                             vecBaselineError)
        [checkMeanEff, _unused] = _lrp_meanPcrEff(tarGroup, vecTarget, checkPcrEff,
                                                  vecSkipSample, vecNoPlateau, vecShortLogLin)
        if checkMeanEff < 1.001:
            skipGroup = True

    if skipGroup:
        if tarGroup is None:
            threshold[0] = (0.5 * np.round(1000 * upWin[0]) / 1000)
        else:
            threshold[tarGroup] = (0.5 * np.round(1000 * upWin[tarGroup]) / 1000)

    if not skipGroup:
        foldWidth = np.log10(np.power(checkMeanEff, pointsInWoL))
        counter = -1
        maxVarEff = 0.0
        maxVarEffStep = -1
        lastUpWin = 2 + maxFluorWin
        while True:
            counter += 1
            step = np.log10(checkMeanEff)
            newUpWin = maxFluorWin - counter * stepSize * step
            if newUpWin < lastUpWin:
                upWin, lowWin = _lrp_setLogWin(tarGroup, newUpWin, foldWidth, upWin, lowWin, maxFluorTotal, minFluorTotal)
                _unused, _unused2, checkPcrEff, _unused3, _unused4, _unused5 = _lrp_allParamInWindow(fluor, tarGroup,
                                                                                                     vecTarget, indMeanX,
                                                                                                     indMeanY, pcrEff,
                                                                                                     nNulls, nInclu,
                                                                                                     correl,
                                                                                                     startCyc, startCycFix, stopCyc,
                                                                                                     upWin, lowWin,
                                                                                                     vecNoAmplification,
                                                                                                     vecBaselineError)
                [checkMeanEff, _unused] = _lrp_meanPcrEff(tarGroup, vecTarget, checkPcrEff,
                                                          vecSkipSample, vecNoPlateau, vecShortLogLin)
                foldWidth = np.log10(np.power(checkMeanEff, pointsInWoL))
                if foldWidth < 0.5:
                    foldWidth = 0.5  # to avoid width = 0 above stopCyc
                upWin, lowWin = _lrp_setLogWin(tarGroup, newUpWin, foldWidth, upWin, lowWin, maxFluorTotal, minFluorTotal)
                _unused, _unused2, checkPcrEff, _unused3, _unused4, _unused5 = _lrp_allParamInWindow(fluor, tarGroup,
                                                                                                     vecTarget, indMeanX,
                                                                                                     indMeanY, pcrEff,
                                                                                                     nNulls, nInclu,
                                                                                                     correl,
                                                                                                     startCyc, startCycFix, stopCyc,
                                                                                                     upWin, lowWin,
                                                                                                     vecNoAmplification,
                                                                                                     vecBaselineError)
                [checkMeanEff, checkVarEff] = _lrp_meanPcrEff(tarGroup, vecTarget, checkPcrEff,
                                                              vecSkipSample, vecNoPlateau, vecShortLogLin)
                if checkVarEff > 0.0:
                    memVarEff[counter] = np.sqrt(checkVarEff) / checkMeanEff
                else:
                    memVarEff[counter] = 0.0
                if checkVarEff > maxVarEff:
                    maxVarEff = checkVarEff
                    maxVarEffStep = counter
                memUpWin[counter] = newUpWin
                memFoldWidth[counter] = foldWidth
                lastUpWin = newUpWin
            else:
                checkVarEff = 0.0

            if counter >= 60 or newUpWin - foldWidth / (pointsInWoL / 2.0) < minFluorLim or checkVarEff < 0.00000000001:
                break

        # corrections: start
        if checkVarEff < 0.00000000001:
            counter -= 1  # remove window with vareff was 0.0

        validSteps = -1
        while True:
            validSteps += 1
            if memVarEff[validSteps] < 0.000001:
                break
        validSteps -= 1  # i = number of valid steps

        minSmooth = memVarEff[0]
        minStep = 0  # default top window

        # next 3 if conditions on i: added to correct smoothing
        if validSteps == 0:
            minStep = 0

        if 0 < validSteps < 4:
            n = -1
            while True:
                n += 1
                if memVarEff[n] < minSmooth:
                    minSmooth = memVarEff[n]
                    minStep = n
                if n == validSteps:
                    break
        if validSteps >= 4:
            n = 0
            while True:
                n += 1
                smoothVar = 0.0
                for m in range(n - 1, n + 2):
                    smoothVar = smoothVar + memVarEff[m]
                smoothVar = smoothVar / 3.0
                if smoothVar < minSmooth:
                    minSmooth = smoothVar
                    minStep = n

                if n >= validSteps - 1 or n > maxVarEffStep:
                    break
        # corrections: stop

        # Calculate the final values again
        upWin, lowWin = _lrp_setLogWin(tarGroup, memUpWin[minStep], memFoldWidth[minStep],
                                       upWin, lowWin, maxFluorTotal, minFluorTotal)
        if tarGroup is None:
            threshold[0] = (0.5 * np.round(1000 * upWin[0]) / 1000)
        else:
            threshold[tarGroup] = (0.5 * np.round(1000 * upWin[tarGroup]) / 1000)

        indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl = _lrp_allParamInWindow(fluor, tarGroup, vecTarget,
                                                                                   indMeanX, indMeanY, pcrEff, nNulls,
                                                                                   nInclu, correl, startCyc, startCycFix, 
                                                                                   stopCyc, upWin, lowWin,
                                                                                   vecNoAmplification, vecBaselineError)
        for aRow in range(0, len(pcrEff)):
            if tarGroup is None or tarGroup == vecTarget[aRow]:
                if (not (vecSkipSample[aRow] or vecNoPlateau[aRow] or vecShortLogLin[aRow])) and pcrEff[aRow] > 1.0:
                    vecIsUsedInWoL[aRow] = True
                else:
                    vecIsUsedInWoL[aRow] = False

    return indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL


def _lrp_assignNoPlateau(fluor, tarGroup, vecTarget, pointsInWoL, indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl,
                         upWin, lowWin, maxFluorTotal, minFluorTotal, stopCyc, startCyc, startCycFix, threshold,
                         vecNoAmplification, vecBaselineError, vecSkipSample, vecNoPlateau, vecShortLogLin, vecIsUsedInWoL):
    """Assign no plateau again and possibly recalculate WoL if new no plateau was found.

    Args:
        fluor: The array with the fluorescence values
        tarGroup: The target number
        vecTarget: The vector with the targets numbers
        pointsInWoL: The number of points in the window
        indMeanX: The vector with the x mean position
        indMeanY: The vector with the y mean position
        pcrEff: The array with the PCR efficiencies
        nNulls: The array with the calculated nNulls
        nInclu: The array with the calculated nInclu
        correl: The array with the calculated correl
        upWin: The upper limit of the window
        lowWin: The lower limit of the window
        maxFluorTotal: The maximum fluorescence over all rows
        minFluorTotal: The minimum fluorescence over all rows
        stopCyc: The vector with the stop cycle of the log lin phase
        startCyc: The vector with the start cycle of the log lin phase
        threshold: The threshold fluorescence
        vecNoAmplification: True if there is a amplification error
        vecBaselineError: True if there is a baseline error
        vecSkipSample: Skip the sample
        vecNoPlateau: True if there is no plateau
        vecShortLogLin: True indicates a short log lin phase
        vecIsUsedInWoL: True if used in the WoL

    Returns:
        The values indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL, vecNoPlateau.
    """
    newNoPlateau = False
    for aRow in range(0, fluor.shape[0]):
        if (tarGroup is None or tarGroup == vecTarget[aRow]) and not (vecNoAmplification[aRow] or
                                                                      vecBaselineError[aRow] or
                                                                      vecNoPlateau[aRow]):
            expectedFluor = nNulls[aRow] * np.power(pcrEff[aRow], fluor.shape[1])
            if expectedFluor / fluor[aRow, fluor.shape[1] - 1] < 5:
                newNoPlateau = True
                vecNoPlateau[aRow] = True

    if newNoPlateau:
        indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL = _lrp_setWoL(fluor, tarGroup, vecTarget,
                                                                                                                   pointsInWoL, indMeanX, indMeanY, pcrEff,
                                                                                                                   nNulls, nInclu, correl, upWin,
                                                                                                                   lowWin, maxFluorTotal, minFluorTotal,
                                                                                                                   stopCyc, startCyc, startCycFix, threshold,
                                                                                                                   vecNoAmplification,
                                                                                                                   vecBaselineError,
                                                                                                                   vecSkipSample, vecNoPlateau,
                                                                                                                   vecShortLogLin, vecIsUsedInWoL)

    return indMeanX, indMeanY, pcrEff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL, vecNoPlateau


def _lrp_removeOutlier(data, vecNoPlateau, alpha=0.05):
    """A function which calculates the skewness and Grubbs test to identify outliers ignoring nan.

    Args:
        data: The numpy array with the data
        vecNoPlateau: The vector of samples without plateau.
        alpha: The the significance level (default 0.05)

    Returns:
        The a bool array with the removed outliers set true.
    """

    oData = np.copy(data)
    oLogic = np.zeros(data.shape, dtype=np.bool_)
    loopOn = True

    while loopOn:
        count = np.count_nonzero(~np.isnan(oData))
        if count < 3:
            loopOn = False
        else:
            mean = np.nanmean(oData)
            std = np.nanstd(oData, ddof=1)
            skewness = scp.skew(oData, bias=False, nan_policy='omit')
            skewness_SE = np.sqrt((6 * count * (count - 1)) / ((count - 2) * (count + 1) * (count + 3)))
            skewness_t = np.abs(skewness) / skewness_SE
            skewness_P = scp.t.sf(skewness_t, df=np.power(10, 10)) * 2

            if skewness_P < alpha / 2.0:
                # It's skewed!
                grubbs_t = scp.t.ppf(1 - (alpha / count) / 2, (count - 2))
                grubbs_Gcrit = ((count - 1) / np.sqrt(count)) * np.sqrt(np.power(grubbs_t, 2) /
                                                                        ((count - 2) + np.power(grubbs_t, 2)))
                if skewness > 0.0:
                    data_max = np.nanmax(oData)
                    grubbs_res = (data_max - mean) / std
                    max_pos = np.nanargmax(oData)
                    if grubbs_res > grubbs_Gcrit:
                        # It's a true outlier
                        oData[max_pos] = np.nan
                        oLogic[max_pos] = True
                    else:
                        if vecNoPlateau[max_pos]:
                            # It has no plateau
                            oData[max_pos] = np.nan
                            oLogic[max_pos] = True
                        else:
                            loopOn = False
                else:
                    data_min = np.nanmin(oData)
                    grubbs_res = (mean - data_min) / std
                    min_pos = np.nanargmin(oData)
                    if grubbs_res > grubbs_Gcrit:
                        # It's a true outlier
                        oData[min_pos] = np.nan
                        oLogic[min_pos] = True
                    else:
                        if vecNoPlateau[min_pos]:
                            # It has no plateau
                            oData[min_pos] = np.nan
                            oLogic[min_pos] = True
                        else:
                            loopOn = False
            else:
                loopOn = False
    return oLogic


def _mca_smooth(tempList, rawFluor):
    """A function to smooth the melt curve date based on Friedmans supersmoother.
       # https://www.slac.stanford.edu/pubs/slacpubs/3250/slac-pub-3477.pdf

    Args:
        tempList:
        rawFluor: The numpy array with the raw data

    Returns:
        The numpy array with the smoothed data.
    """

    span_s = 0.05
    span_m = 0.2
    span_l = 0.5

    smoothFluor = np.zeros(rawFluor.shape, dtype=np.float64)

    padTemp = np.append(0.0, tempList)

    zeroPad = np.zeros((rawFluor.shape[0], 1), dtype=np.float64)
    padFluor = np.append(zeroPad, rawFluor, axis=1)
    n = len(padTemp) - 1

    # Find the increase in x from 0.25 to 0.75 over the total range
    firstQuarter = int(0.5 + n / 4)
    thirdQuarter = 3 * firstQuarter
    scale = -1.0
    while scale <= 0.0:
        if thirdQuarter < n:
            thirdQuarter += 1
        if firstQuarter > 1:
            firstQuarter -= 1
        scale = padTemp[thirdQuarter] - padTemp[firstQuarter]
    vsmlsq = 0.0001 * scale * 0.0001 * scale

    countUp = 0
    for fluor in padFluor:
        [res_s_a, res_s_t] = _mca_sub_smooth(padTemp, fluor, span_s, vsmlsq, True)
        [res_s_b, _unused] = _mca_sub_smooth(padTemp, res_s_t, span_m, vsmlsq, False)
        [res_s_c, res_s_t] = _mca_sub_smooth(padTemp, fluor, span_m, vsmlsq, True)
        [res_s_d, _unused] = _mca_sub_smooth(padTemp, res_s_t, span_m, vsmlsq, False)
        [res_s_e, res_s_t] = _mca_sub_smooth(padTemp, fluor, span_l, vsmlsq, True)
        [res_s_f, _unused] = _mca_sub_smooth(padTemp, res_s_t, span_m, vsmlsq, False)

        res_s_fin = np.zeros(res_s_a.shape, dtype=np.float64)
        for thirdQuarter in range(1, n + 1):
            resmin = 1.0e20
            if res_s_b[thirdQuarter] < resmin:
                resmin = res_s_b[thirdQuarter]
                res_s_fin[thirdQuarter] = span_s
            if res_s_d[thirdQuarter] < resmin:
                resmin = res_s_d[thirdQuarter]
                res_s_fin[thirdQuarter] = span_m
            if res_s_f[thirdQuarter] < resmin:
                res_s_fin[thirdQuarter] = span_l

        [res_s_bb, _unused] = _mca_sub_smooth(padTemp, res_s_fin, span_m, vsmlsq, False)

        res_s_cc = np.zeros(res_s_a.shape, dtype=np.float64)
        for thirdQuarter in range(1, n + 1):
            # compare res_s_bb with spans[] and make sure the no res_s_bb[] is below span_s or above span_l
            if res_s_bb[thirdQuarter] <= span_s:
                res_s_bb[thirdQuarter] = span_s
            if res_s_bb[thirdQuarter] >= span_l:
                res_s_bb[thirdQuarter] = span_l
            f = res_s_bb[thirdQuarter] - span_m
            if f >= 0.0:
                # in case res_s_bb[] is higher than span_m: calculate res_s_cc[] from res_s_c and res_s_e
                # using linear interpolation between span_l and span_m
                f = f / (span_l - span_m)
                res_s_cc[thirdQuarter] = (1.0 - f) * res_s_c[thirdQuarter] + f * res_s_e[thirdQuarter]
            else:
                # in case res_s_bb[] is less than span_m: calculate res_s_cc[] from res_s_c and res_s_a
                # using linear interpolation between span_s and span_m
                f = -f / (span_m - span_s)
                res_s_cc[thirdQuarter] = (1.0 - f) * res_s_c[thirdQuarter] + f * res_s_a[thirdQuarter]

        # final smoothing of combined optimally smoothed values in res_s_cc[] into smo[]
        [res_s_t, _unused] = _mca_sub_smooth(padTemp, res_s_cc, span_s, vsmlsq, False)
        smoothFluor[countUp] = res_s_t[1:]
        countUp += 1

    return smoothFluor


def _mca_sub_smooth(temperature, fluor, span, vsmlsq, saveVarianceData):
    """A function to smooth the melt curve date based on Friedmans supersmoother.
       # https://www.slac.stanford.edu/pubs/slacpubs/3250/slac-pub-3477.pdf

    Args:
        temperature:
        fluor: The numpy array with the raw data
        span: The selected span
        vsmlsq: The width
        saveVarianceData: Sava variance data

    Returns:
        [smoothData[], varianceData[]]  where smoothData[] contains smoothed data,
        varianceData[] contains residuals scaled to variance.
    """

    n = len(temperature) - 1
    smoothData = np.zeros(len(temperature), dtype=np.float64)
    varianceData = np.zeros(len(temperature), dtype=np.float64)

    windowSize = int(0.5 * span * n + 0.6)
    if windowSize < 2:
        windowSize = 2
    windowStop = 2 * windowSize + 1  # range of smoothing window

    xm = temperature[1]
    ym = fluor[1]
    tempVar = 0.0
    fluorVar = 0.0

    for i in range(2, windowStop + 1):
        xm = ((i - 1) * xm + temperature[i]) / i
        ym = ((i - 1) * ym + fluor[i]) / i
        tmp = i * (temperature[i] - xm) / (i - 1)
        tempVar += tmp * (temperature[i] - xm)
        fluorVar += tmp * (fluor[i] - ym)

    fbw = windowStop
    for j in range(1, n + 1):  # Loop through all
        windowStart = j - windowSize - 1
        windowEnd = j + windowSize
        if not (windowStart < 1 or windowEnd > n):
            tempStart = temperature[windowStart]
            tempEnd = temperature[windowEnd]
            fbo = fbw
            fbw = fbw - 1.0
            tmp = 0.0
            if fbw > 0.0:
                xm = (fbo * xm - tempStart) / fbw
            if fbw > 0.0:
                ym = (fbo * ym - fluor[windowStart]) / fbw
            if fbw > 0.0:
                tmp = fbo * (tempStart - xm) / fbw
            tempVar = tempVar - tmp * (tempStart - xm)
            fluorVar = fluorVar - tmp * (fluor[windowStart] - ym)

            fbo = fbw
            fbw = fbw + 1.0
            tmp = 0.0
            if fbw > 0.0:
                xm = (fbo * xm + tempEnd) / fbw
            if fbw > 0.0:
                ym = (fbo * ym + fluor[windowEnd]) / fbw
            if fbo > 0.0:
                tmp = fbw * (tempEnd - xm) / fbo
            tempVar = tempVar + tmp * (tempEnd - xm)
            fluorVar = fluorVar + tmp * (fluor[windowEnd] - ym)

        if tempVar > vsmlsq:
            smoothData[j] = (temperature[j] - xm) * fluorVar / tempVar + ym  # contains smoothed data
        else:
            smoothData[j] = ym  # contains smoothed data

        if saveVarianceData:
            h = 0.0
            if fbw > 0.0:
                h = 1.0 / fbw
            if tempVar > vsmlsq:
                h = h + (temperature[j] - xm) * (temperature[j] - xm) / tempVar

            if 1.0 - h > 0.0:
                varianceData[j] = abs(fluor[j] - smoothData[j]) / (1.0 - h)  # contains residuals scaled to variance
            else:
                if j > 1:
                    varianceData[j] = varianceData[j - 1]  # contains residuals scaled to variance
                else:
                    varianceData[j] = 0.0

    return [smoothData, varianceData]


def _mca_linReg(xIn, yUse, start, stop):
    """A function which calculates the slope or the intercept by linear regression.

    Args:
        xIn: The numpy array of the temperatures
        yUse: The numpy array that contains the fluorescence

    Returns:
        An array with the slope and intercept.
    """

    counts = np.ones(yUse.shape)
    xUse = xIn.copy()
    xUse[np.isnan(yUse)] = 0
    counts[np.isnan(yUse)] = 0
    myStop = stop + 1

    tempSqared = xUse * xUse
    tempFluor = xUse * yUse

    sumCyc = np.nansum(xUse[:, start:myStop], axis=1)
    sumFluor = np.nansum(yUse[:, start:myStop], axis=1)
    sumCycSquared = np.nansum(tempSqared[:, start:myStop], axis=1)
    sumCycFluor = np.nansum(tempFluor[:, start:myStop], axis=1)
    n = np.nansum(counts[:, start:myStop], axis=1)

    ssx = sumCycSquared - (sumCyc * sumCyc) / n
    sxy = sumCycFluor - (sumCyc * sumFluor) / n

    slope = sxy / ssx
    intercept = (sumFluor / n) - slope * (sumCyc / n)
    return [slope, intercept]


def _pco_fixPlateMatix(mat, row, col):
    """A function which calculates the missing values for the given position.

    Args:
        mat: The matrix to correct
        row: The row position of the value to correct
        col: The column position of the value to correct

    Returns:
        None, changes the values in the matrix.
    """

    matSize = len(mat)

    # Return if there are useful no values
    maxVal = -10
    for curr in range(0, matSize):
        maxVal = max(maxVal, mat[row][curr])
        maxVal = max(maxVal, mat[curr][col])
    if maxVal <= 0.0:
        return

    # Calc the missing values
    finSum = 0.0
    finNum = 0
    for mCol in range(0, matSize):
        if mCol == col:
            continue
        colSum = 0.0
        colNum = 0
        for mRow in range(0, matSize):
            if mRow == row:
                continue
            if mat[mRow][mCol] <= 0.0:
                continue
            if mat[mRow][col] <= 0.0:
                continue
            colSum += math.log(mat[mRow][mCol] / mat[mRow][col])
            colNum += 1
        if colNum > 0:
            colGeo = math.exp(colSum / colNum)
            if colGeo > 0:
                finSum += math.log(mat[row][mCol] / colGeo)
                finNum += 1
    finFact = -1.0
    if finNum > 0:
        finFact = math.exp(finSum / finNum)
    if finNum > 0:
        mat[row][col] = finFact
        mat[col][row] = 1.0 / finFact


def _cleanErrorString(inStr, cleanStyle):
    outStr = ";"
    inStr += ";"
    if cleanStyle == "melt":
        outStr = inStr.replace('several products with different melting temperatures detected', '')
        outStr = outStr.replace('product with different melting temperatures detected', '')
        outStr = outStr.replace('no product with expected melting temperature', '')
        outStr = outStr.replace('product detected in negative control', '')
    else:
        strList = inStr.split(";")
        knownWarn = ["amplification in negative control", "plateau in negative control",
                     "no amplification in positive control", "baseline error in positive control",
                     "instable baseline", "instable baseline in positive control",
                     "no plateau in positive control", "noisy sample in positive control",
                     "Ncopy > 100 Mio", "Ncopy < 10", "no indiv PCR eff can be calculated",
                     "PCR efficiency outlier", "no amplification", "baseline error", "no plateau",
                     "noisy sample", "Cq < 10", "N0 unreliable", "Cq > 34", "Cq too high"]
        for ele in strList:
            if ele in knownWarn:
                continue
            if re.search(r"^only \d+ values in log phase", ele):
                continue
            if re.search(r"^indiv PCR eff is .+", ele):
                continue
            outStr += ele + ";"

     #   if inStr.find('several products with different melting temperatures detected') >= 0:
     #       outStr += ';several products with different melting temperatures detected;'
     #   if inStr.find('product with different melting temperatures detected') >= 0:
     #       outStr += ';product with different melting temperatures detected;'
     #   if inStr.find('no product with expected melting temperature') >= 0:
     #       outStr += ';no product with expected melting temperature;'

    outStr = re.sub(r';+', ';', outStr)
    return outStr


def _numpyTwoAxisSave(var, fileName):
    with np.printoptions(precision=3, suppress=True):
        np.savetxt(fileName, var, fmt='%.6f', delimiter='\t', newline='\n')


def _getXMLDataType():
    return ["tar", "cq", "N0", "Ncopy", "ampEffMet", "ampEff", "ampEffSE", "corrF",
            "corrP", "corrCq", "meltTemp", "excl", "note", "adp", "mdp", "endPt",
            "bgFluor", "quantFluor"]


def _getXMLPartitionDataType():
    return ["tar", "excluded", "note", "pos", "neg", "undef", "excl", "conc"]

def runStatistics(statTarGroup, parametric, translateGrp, statAlpha=0.05):
    try:
        statAlpha = float(statAlpha)
    except ValueError:
        raise RdmlError('Error: Alpha must be a number.')
    if statAlpha <= 0.0:
        raise RdmlError('Error: Alpha must be > 0.0.')
    if statAlpha >= 1.0:
        raise RdmlError('Error: Alpha must be < 1.0.')
    sym = '#'
    if statAlpha <= 0.05:
        sym = '*'
    if statAlpha <= 0.01:
        sym = '**'
    if statAlpha <= 0.001:
        sym = '***'

    ret = {}
    ret["multi comparison"] = ""
    ret["alpha"] = statAlpha
    ret["df"] = len(statTarGroup) - 1
    if len(statTarGroup) < 2:
        ret["test name"] = "no test possible"
        ret["stat name"] = "only 1 group"
        ret["stat val"] = -1.0
        ret["p val"] = -1.0
        return ret

    if _string_to_bool(parametric, triple=False):
        if len(statTarGroup) == 2:
            ret["test name"] = "T-test"
            ret["stat name"] = "t-statistic"
            ret["df"] = len(statTarGroup[0]) + len(statTarGroup[1]) - 2
            ret["stat val"], ret["p val"] = scp.ttest_ind(statTarGroup[0], statTarGroup[1])
        else:
            ret["test name"] = "One-way ANOVA"
            ret["stat name"] = "F statistic"
            ret["stat val"], ret["p val"] = scp.f_oneway(*statTarGroup)

            if ret["p val"] > statAlpha:
                return ret
            gridGroup = []
            gridValue = []
            grid_Y2 = []
            gridGroup = []
            groupRes = []
            groupIncreasMean = []
            grpGridRanksum = []
            statres_total_n = 0
            grpGridObsDiff = []
            grpGridGrpRange = []
            grpGridCriticQ = []
            grpGridCriticDiff = []
            grpGridMultiComp = []
            # Prepare Data
            for row in range(0, len(statTarGroup)):
                groupRes.append({})
                groupRes[row]["n"] = 0
                groupRes[row]["median"] = np.median(statTarGroup[row])
                groupRes[row]["sumY"] = 0.0
                groupRes[row]["sumY2"] = 0.0
                groupRes[row]["mean"] = 0.0
                groupRes[row]["usedMean"] = False
                groupRes[row]["var"] = 0.0
                groupRes[row]["SD"] = 0.0
                for col in range(0, len(statTarGroup[row])):
                    gridGroup.append(row)
                    gridValue.append(statTarGroup[row][col])
                    grid_Y2.append(statTarGroup[row][col] * statTarGroup[row][col])
                    groupRes[row]["n"] += 1
                    groupRes[row]["sumY"] += statTarGroup[row][col]
                    groupRes[row]["sumY2"] += statTarGroup[row][col] * statTarGroup[row][col]
            statres_total_n = len(gridValue)
            # Calculate the group results
            stat_quant1 = 0.0
            stat_quant2 = 0.0
            stat_quant3 = 0.0
            stat_quant4 = 0.0
            for row in range(0, len(groupRes)):
                groupRes[row]["sumY2_n"] = groupRes[row]["sumY"] * groupRes[row]["sumY"] / groupRes[row]["n"]
                stat_quant1 += groupRes[row]["sumY"]
                stat_quant2 += groupRes[row]["sumY2"]
                stat_quant3 += groupRes[row]["sumY2_n"]
                groupRes[row]["mean"] = groupRes[row]["sumY"] / groupRes[row]["n"]
                groupRes[row]["var"] = (groupRes[row]["sumY2"] - groupRes[row]["sumY2_n"]) / (groupRes[row]["n"] - 1)
                groupRes[row]["SD"] = np.sqrt(groupRes[row]["var"])
                groupRes[row]["SEM"] = np.sqrt(groupRes[row]["var"]) / np.sqrt(groupRes[row]["n"])
            stat_quant4 = stat_quant1 * stat_quant1 / statres_total_n
            stat_SStotal = stat_quant2 - stat_quant4
            stat_SSgroups = stat_quant3 - stat_quant4
            stat_SSwithin = stat_SStotal - stat_SSgroups

            stat_df_total = statres_total_n - 1
            stat_df_groups = len(statTarGroup) - 1
            stat_df_within = statres_total_n - len(statTarGroup)

            stat_ms_groups = stat_SSgroups / stat_df_groups
            stat_ms_within = stat_SSwithin / stat_df_within

            stat_F = stat_ms_groups / stat_ms_within
            stat_P = 1 - scp.f.cdf(stat_F, stat_df_groups, stat_df_within)
            stat_grp_range = scp.t.ppf(1 - statAlpha / 2.0, stat_df_within)

            # Student-Newman-Keuls multiple comparison of groups
            # Sort by increasing mean
            sortCurrNr = -1
            sortMinVal = 999999999999999
            sortOn = True
            while sortOn:
                for row in range(0, len(groupRes)):
                    if not groupRes[row]["usedMean"]:
                        if groupRes[row]["mean"] < sortMinVal:
                            sortCurrNr = row
                            sortMinVal = groupRes[row]["mean"]
                if sortMinVal != 999999999999999:
                    groupIncreasMean.append(sortCurrNr)
                    groupRes[sortCurrNr]["usedMean"] = True
                    sortMinVal = 999999999999999
                else:
                    sortOn = False
            # Calculate observed difference
            for row in range(0, len(groupIncreasMean)):
                grpGridObsDiff.append([])
                grpGridGrpRange.append([])
                for col in range(0, len(groupIncreasMean)):
                    grpGridObsDiff[row].append(np.abs(groupRes[groupIncreasMean[row]]["mean"] - groupRes[groupIncreasMean[col]]["mean"]))
                    grpGridGrpRange[row].append(np.abs(row - col) + 1)
            # Calculate critical Q from studentized range distribution
            for row in range(0, len(groupIncreasMean)):
                grpGridCriticQ.append([])
                for col in range(0, len(groupIncreasMean)):
                    if grpGridGrpRange[row][col] == 1:
                        grpGridCriticQ[row].append("")
                    elif grpGridGrpRange[row][col] > stat_df_within:
                        grpGridCriticQ[row].append("")
                    else:
                        grpGridCriticQ[row].append((stat_grp_range * (np.log(grpGridGrpRange[row][col] - 1) *
                                                                      (0.8843 - 0.2368 * stat_grp_range -
                                                                       1.214 / stat_df_within +
                                                                       1.208 * stat_grp_range / stat_df_within) +
                                                                       np.sqrt(2))))
            # Calculate critical difference
            for row in range(0, len(statTarGroup)):
                grpGridCriticDiff.append([])
                for col in range(0, len(statTarGroup)):
                    if grpGridCriticQ[row][col] == "":
                        grpGridCriticDiff[row].append("")
                    else:
                        grpGridCriticDiff[row].append(grpGridCriticQ[row][col] *
                                                      np.sqrt(stat_ms_within) *
                                                      np.sqrt((groupRes[groupIncreasMean[row]]["n"] + groupRes[groupIncreasMean[col]]["n"]) /
                                                                (2.0 * groupRes[groupIncreasMean[row]]["n"] * groupRes[groupIncreasMean[col]]["n"])))
            # Calculate multiple comparison result
            grpGridMultiComp.append(["group", ""])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp[0].append(translateGrp[groupIncreasMean[row]])
            grpGridMultiComp.append(["", "mean"])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp[1].append(groupRes[groupIncreasMean[row]]["mean"])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp.append([translateGrp[groupIncreasMean[row]], groupRes[groupIncreasMean[row]]["mean"]])
                for col in range(0, len(statTarGroup)):
                    if grpGridCriticDiff[row][col] == "":
                        grpGridMultiComp[row + 2].append("-")
                    elif grpGridObsDiff[row][col] > grpGridCriticDiff[row][col]:
                        grpGridMultiComp[row + 2].append(sym)
                    else:
                        grpGridMultiComp[row + 2].append("ns")
            ret["multi comparison"] = grpGridMultiComp
    else:
        if len(statTarGroup) == 2:
            ret["test name"] = "Mann-Whitney U rank test"
            ret["stat name"] = "U statistic"
            ret["stat val"], ret["p val"] = scp.mannwhitneyu(statTarGroup[0], statTarGroup[1], alternative='two-sided', use_continuity=False, method='asymptotic')
        else:
            ret["test name"] = "Kruskal-Wallis H-test"
            ret["stat name"] = "H statistic"
            ret["stat val"], ret["p val"] = scp.kruskal(*statTarGroup)

            if ret["p val"] > statAlpha:
                return ret
            gridGroup = []
            gridValue = []
            scoreGrid = []
            scoreRank = []
            groupRes = []
            grpGridRanksum = []
            grpGridCriticDiff = []
            grpGridMultiComp = []
            statres_total_n = 0
            statres_s1 = 0.0
            statres_sqrsum = 0.0
            # Prepare Data
            for row in range(0, len(statTarGroup)):
                groupRes.append({})
                groupRes[row]["n"] = 0
                groupRes[row]["median"] = np.median(statTarGroup[row])
                groupRes[row]["rsum"] = 0.0
                groupRes[row]["rsum2_n"] = 0.0
                for col in range(0, len(statTarGroup[row])):
                    gridGroup.append(row)
                    gridValue.append(statTarGroup[row][col])
            statres_total_n = len(gridValue)
            # Create the score grid filled
            for row in range(0, len(gridValue)):
                scoreGrid.append([])
                for col in range(0, len(gridValue)):
                    if gridValue[row] == gridValue[col]:
                        scoreGrid[row].append(0.5)
                    elif gridValue[row] < gridValue[col]:
                        scoreGrid[row].append(0.0)
                    elif gridValue[row] > gridValue[col]:
                        scoreGrid[row].append(1.0)
            # Calculate the sums
            for row in range(0, len(scoreGrid)):
                sumUp = 0.5
                for col in range(0, len(scoreGrid[row])):
                    sumUp += scoreGrid[row][col]
                scoreRank.append(sumUp)
                groupRes[gridGroup[row]]["n"] += 1
                groupRes[gridGroup[row]]["rsum"] += sumUp
                statres_sqrsum += sumUp * sumUp
            # Calculate the group results
            for row in range(0, len(groupRes)):
                groupRes[row]["rsum2_n"] = groupRes[row]["rsum"] * groupRes[row]["rsum"] / groupRes[row]["n"]
                statres_s1 += groupRes[row]["rsum2_n"]
            # Final calculations
            statres_n3 = statres_total_n * (statres_total_n + 1) * (statres_total_n + 1) / 4
            statres_s2 = (1 / (statres_total_n - 1 )) * (statres_sqrsum - statres_n3)
            statres_tt = (1 / statres_s2) * (statres_s1 - statres_n3)
            statres_df = len(statTarGroup) - 1
            statres_prob = 1 - scp.chi2.cdf(statres_tt, df=statres_df)
            statres_t = scp.t.ppf(1 - statAlpha / 2.0, statres_total_n - len(statTarGroup))
            statres_SD = statres_t * np.sqrt(( statres_s2 * (( statres_total_n - 1 - statres_tt) / (statres_total_n - len(statTarGroup)))))
            # Calculate mean ranksum difference
            for row in range(0, len(statTarGroup)):
                grpGridRanksum.append([])
                for col in range(0, len(statTarGroup)):
                    grpGridRanksum[row].append(np.abs(groupRes[row]["rsum"] / groupRes[row]["n"] - groupRes[col]["rsum"] / groupRes[col]["n"]))
            # Calculate critical difference
            for row in range(0, len(statTarGroup)):
                grpGridCriticDiff.append([])
                for col in range(0, len(statTarGroup)):
                    if grpGridRanksum[row][col] <= 0.0:
                        grpGridCriticDiff[row].append("")
                    else:
                        grpGridCriticDiff[row].append(statres_SD * np.sqrt(1 / groupRes[row]["n"] + 1 / groupRes[col]["n"]))
            # Calculate multiple comparison result
            grpGridMultiComp.append(["group", ""])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp[0].append(translateGrp[row])
            grpGridMultiComp.append(["", "median"])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp[1].append(groupRes[row]["median"])
            for row in range(0, len(statTarGroup)):
                grpGridMultiComp.append([translateGrp[row], groupRes[row]["median"]])
                for col in range(0, len(statTarGroup)):
                    if grpGridCriticDiff[row][col] == "":
                        grpGridMultiComp[row + 2].append("-")
                    elif grpGridRanksum[row][col] > grpGridCriticDiff[row][col]:
                        grpGridMultiComp[row + 2].append(sym)
                    else:
                        grpGridMultiComp[row + 2].append("ns")
            ret["multi comparison"] = grpGridMultiComp
    return ret

def webAppRunStatistics(data, parametric=False, statAlpha=0.05, seperator='\t', replaceComma=True):
    replaceComma = _string_to_bool(replaceComma, triple=False)
    goodData = []
    foundGrp = {}
    translateGrp = {}
    posGrp = 0

    dff = io.StringIO(data)
    rawData = list(csv.reader(dff, delimiter=seperator))
    for row in range(0, len(rawData)):
        if len(rawData[row]) < 2:
            continue
        if rawData[row][0] not in foundGrp:
            foundGrp[rawData[row][0]] = posGrp
            translateGrp[posGrp] = rawData[row][0]
            goodData.append([])
            posGrp += 1
        for col in range(1, len(rawData[row])):
            cell = rawData[row][col]
            if replaceComma:
                cell = re.sub(r"\.", "", cell)
                cell = re.sub(r",", r"\.", cell)
            goodData[foundGrp[rawData[row][0]]].append(float(cell))

    return runStatistics(goodData, parametric, translateGrp, statAlpha)

def standardCurveStats(data, noCq=False):
    # The data are as data{target}{"raw"}{"standard conc as string"}["float, "float"]
    targets = sorted(list(data.keys()))

    for tar in range(0, len(targets)):
        # Values indicating failure
        data[targets[tar]]["cq_eff"] = -1.0
        data[targets[tar]]["mean_min"] = -1.0
        data[targets[tar]]["mean_max"] = -1.0
        data[targets[tar]]["bias_as_ratio"] = -1.0
        data[targets[tar]]["expected_min"] = -1.0
        data[targets[tar]]["expected_max"] = -1.0
        data[targets[tar]]["expected_bias_as_ratio"] = -1.0
        data[targets[tar]]["slope_bias"] = -1.0
        data[targets[tar]]["correlation_R"] = -1.0
        data[targets[tar]]["MSderivation"] = -1.0  # Linearity
        data[targets[tar]]["MSwithin"] = -1.0  # Precision (Reproducibility)
        data[targets[tar]]["cq_ms_within_linregpcr"] = -1.0  # MS within LinRegPCR
        data[targets[tar]]["cq_ms_ratio"] = -1.0  # Ratio
        data[targets[tar]]["dif_detectable_diff"] = -1.0  # Resolution (detectable difference)

        # Remove negative and 0 values to be log save
        for stdval in data[targets[tar]]["raw"]:
            for pos in range(0, len(data[targets[tar]]["raw"][stdval])):
                if data[targets[tar]]["raw"][stdval][pos] <= 0.0:
                    del data[targets[tar]]["raw"][stdval][pos]

        # Extract the information about the standard
        data[targets[tar]]["std"] = {}
        data[targets[tar]]["std"]["val"] = []
        data[targets[tar]]["std"]["str"] = []
        for stdval in data[targets[tar]]["raw"]:
            stdF = float(re.sub(r"[^0-9\-\.]", "", stdval))
            if len(data[targets[tar]]["std"]["val"]) == 0:
                data[targets[tar]]["std"]["val"].append(stdF)
                data[targets[tar]]["std"]["str"].append(stdval)
            else:
                for pos in range(0, len(data[targets[tar]]["std"]["val"])):
                    if data[targets[tar]]["std"]["val"][pos] > stdF:
                        data[targets[tar]]["std"]["val"].insert(pos, stdF)
                        data[targets[tar]]["std"]["str"].insert(pos, stdval)
                        break
                    if pos == len(data[targets[tar]]["std"]["val"]) - 1:
                        data[targets[tar]]["std"]["val"].append(stdF)
                        data[targets[tar]]["std"]["str"].append(stdval)

        # Ignore targets without <2 values
        if len(data[targets[tar]]["std"]["str"]) < 1:
            continue

        # Calculate Min / Max Vals
        maxPosition = len(data[targets[tar]]["std"]["val"]) - 1
        if data[targets[tar]]["std"]["val"][0] == data[targets[tar]]["std"]["val"][maxPosition]:
            continue
        data[targets[tar]]["var_groups"] = len(data[targets[tar]]["raw"])
        data[targets[tar]]["mean_min"] = np.mean(data[targets[tar]]["raw"][data[targets[tar]]["std"]["str"][0]])
        data[targets[tar]]["mean_max"] = np.mean(data[targets[tar]]["raw"][data[targets[tar]]["std"]["str"][maxPosition]])
        data[targets[tar]]["bias_as_ratio"] = data[targets[tar]]["mean_max"] / data[targets[tar]]["mean_min"]
        data[targets[tar]]["expected_min"] = data[targets[tar]]["std"]["val"][0]
        data[targets[tar]]["expected_max"] = data[targets[tar]]["std"]["val"][maxPosition]
        data[targets[tar]]["expected_bias_as_ratio"] = data[targets[tar]]["expected_max"] / data[targets[tar]]["expected_min"]

       # Ignore targets without <2 values
        if len(data[targets[tar]]["std"]["str"]) < 2:
            continue

        # Calulate Log and scaled data
        data[targets[tar]]["std"]["scaled"] = []
        data[targets[tar]]["std"]["log"] = []
        data[targets[tar]]["std"]["scaled_log"] = []
        for pos in range(0, len(data[targets[tar]]["std"]["val"])):
            data[targets[tar]]["std"]["scaled"].append(data[targets[tar]]["std"]["val"][pos] / data[targets[tar]]["expected_max"])
            data[targets[tar]]["std"]["log"].append(np.log10(data[targets[tar]]["std"]["val"][pos]))
            data[targets[tar]]["std"]["scaled_log"].append(np.log10(data[targets[tar]]["std"]["scaled"][pos]))
        data[targets[tar]]["scaled"] = {}
        data[targets[tar]]["log"] = {}
        data[targets[tar]]["scaled_log"] = {}
        data[targets[tar]]["all_std"] = []
        data[targets[tar]]["all_std_scaled_log"] = []
        data[targets[tar]]["all_data_scaled_log"] = []
        data[targets[tar]]["variance"] = {}
        data[targets[tar]]["mean_log"] = {}
        data[targets[tar]]["SSwithin"] = 0.0
        for concPos in range(0, len(data[targets[tar]]["std"]["val"])):
            concStr = data[targets[tar]]["std"]["str"][concPos]
            concFloat = data[targets[tar]]["std"]["val"][concPos]
            data[targets[tar]]["scaled"][concStr] = []
            data[targets[tar]]["log"][concStr] = []
            data[targets[tar]]["scaled_log"][concStr] = []
            for pos in range(0, len(data[targets[tar]]["raw"][concStr])):
                data[targets[tar]]["scaled"][concStr].append(data[targets[tar]]["raw"][concStr][pos] / data[targets[tar]]["mean_max"])
                data[targets[tar]]["log"][concStr].append(np.log10(data[targets[tar]]["raw"][concStr][pos]))
                data[targets[tar]]["scaled_log"][concStr].append(np.log10(data[targets[tar]]["scaled"][concStr][pos]))

            if len(data[targets[tar]]["scaled_log"][concStr]) > 1:
                data[targets[tar]]["variance"][concStr] = np.var(data[targets[tar]]["scaled_log"][concStr], ddof=1)
                data[targets[tar]]["SSwithin"] += data[targets[tar]]["variance"][concStr] * (len(data[targets[tar]]["scaled_log"][concStr]) - 1.0)
            else:
                data[targets[tar]]["variance"][concStr] = -1.0
            data[targets[tar]]["mean_log"][concStr] = np.mean(data[targets[tar]]["log"][concStr])

            for pos in range(0, len(data[targets[tar]]["scaled_log"][concStr])):
                data[targets[tar]]["all_std"].append(concFloat)
                data[targets[tar]]["all_std_scaled_log"].append(np.log10(concFloat / data[targets[tar]]["expected_max"]))
                data[targets[tar]]["all_data_scaled_log"].append(data[targets[tar]]["scaled_log"][concStr][pos])
        data[targets[tar]]["var_n"] = len(data[targets[tar]]["all_data_scaled_log"])
        data[targets[tar]]["scaled_mean_min"] = np.mean(data[targets[tar]]["scaled"][data[targets[tar]]["std"]["str"][0]])
        data[targets[tar]]["scaled_mean_max"] = np.mean(data[targets[tar]]["scaled"][data[targets[tar]]["std"]["str"][maxPosition]])
        if data[targets[tar]]["var_n"] > 2:
            slope, intercept, r_val, p_val, std_err = scp.linregress(data[targets[tar]]["all_std_scaled_log"], data[targets[tar]]["all_data_scaled_log"])
            data[targets[tar]]["slope_bias"] = slope
            data[targets[tar]]["correlation_R"] = r_val
            data[targets[tar]]["all_variance"] = np.var(data[targets[tar]]["all_data_scaled_log"], ddof=1)
            data[targets[tar]]["SStotal"] = data[targets[tar]]["all_variance"] * (data[targets[tar]]["var_n"] - 1)
            data[targets[tar]]["SSbetween"] = data[targets[tar]]["SStotal"] - data[targets[tar]]["SSwithin"]
            # Calc stexy
            fit = np.polyfit(data[targets[tar]]["all_std_scaled_log"], data[targets[tar]]["all_data_scaled_log"], deg=1)
            y_pred = fit[0] * np.asarray(data[targets[tar]]["all_std_scaled_log"]) + fit[1]
            data[targets[tar]]["SSx"] = np.var(data[targets[tar]]["all_std_scaled_log"], ddof=1) * (data[targets[tar]]["var_n"] - 1.0)
            data[targets[tar]]["stexy"] = (((np.asarray(data[targets[tar]]["all_data_scaled_log"]) - y_pred) ** 2).sum() / (data[targets[tar]]["var_n"] - 2.0)) ** 0.5
            data[targets[tar]]["SSres"] = np.power(data[targets[tar]]["stexy"], 2) * (data[targets[tar]]["var_n"] - 2.0)
            data[targets[tar]]["Sres"] = np.power(data[targets[tar]]["SSres"] / (data[targets[tar]]["var_n"] - 2.0), 0.5)
            data[targets[tar]]["se_of_slope"] = data[targets[tar]]["Sres"] / np.power(data[targets[tar]]["SSx"], 0.5)
            data[targets[tar]]["t"] = np.abs(1 - data[targets[tar]]["slope_bias"]) / data[targets[tar]]["se_of_slope"]
            data[targets[tar]]["p"] = scp.t.sf(data[targets[tar]]["t"], df=(data[targets[tar]]["var_n"] - 2.0)) * 2.0
            data[targets[tar]]["SSregression"] = data[targets[tar]]["SStotal"] - data[targets[tar]]["SSres"]
            data[targets[tar]]["SSderivation"] = data[targets[tar]]["SSres"] - data[targets[tar]]["SSwithin"]
            data[targets[tar]]["MSbetween"] = data[targets[tar]]["SSbetween"] / (data[targets[tar]]["var_groups"] - 1.0)
            data[targets[tar]]["MSregression"] = data[targets[tar]]["SSregression"] / 1
            data[targets[tar]]["MSderivation"] = data[targets[tar]]["SSderivation"] / (data[targets[tar]]["var_groups"] - 2.0)
            data[targets[tar]]["MSwithin"] = data[targets[tar]]["SSwithin"] / (data[targets[tar]]["var_n"] - data[targets[tar]]["var_groups"])
        else:
            data[targets[tar]]["slope_bias"] = -1.0
            data[targets[tar]]["correlation_R"] = -1.0
            data[targets[tar]]["all_variance"] = -1.0
            data[targets[tar]]["SStotal"] = -1.0
            data[targets[tar]]["SSbetween"] = -1.0
            data[targets[tar]]["SSx"] = -1.0
            data[targets[tar]]["stexy"] = -1.0
            data[targets[tar]]["SSres"] = -1.0
            data[targets[tar]]["Sres"] = -1.0
            data[targets[tar]]["se_of_slope"] = -1.0
            data[targets[tar]]["t"] = -1.0
            data[targets[tar]]["p"] = -1.0
            data[targets[tar]]["SSregression"] = -1.0
            data[targets[tar]]["SSderivation"] = -1.0
            data[targets[tar]]["MSbetween"] = -1.0
            data[targets[tar]]["MSregression"] = -1.0
            data[targets[tar]]["MSderivation"] = -1.0
            data[targets[tar]]["MSwithin"] = -1.0

        # Calculate things based on Cq
        if not noCq:
            data[targets[tar]]["all_cq"] = []
            data[targets[tar]]["cq_SSwithin"] = 0.0
            data[targets[tar]]["cq_ss_per_dil"] = {}
            data[targets[tar]]["Cq_F0"] = {}
            cq_check_y = []
            for concPos in range(0, len(data[targets[tar]]["std"]["val"])):
                concStr = data[targets[tar]]["std"]["str"][concPos]
                for pos in range(0, len(data[targets[tar]]["Cq"][concStr])):
                    data[targets[tar]]["all_cq"].append(data[targets[tar]]["Cq"][concStr][pos])
            slope, intercept, r_val, p_val, std_err = scp.linregress(data[targets[tar]]["all_std_scaled_log"], data[targets[tar]]["all_cq"])
            data[targets[tar]]["cq_slope_bias"] = slope
            data[targets[tar]]["cq_eff"] = np.power(10, -1.0 / data[targets[tar]]["cq_slope_bias"])

            for concPos in range(0, len(data[targets[tar]]["std"]["val"])):
                concStr = data[targets[tar]]["std"]["str"][concPos]
                if concStr not in data[targets[tar]]["Cq_F0"]:
                    data[targets[tar]]["Cq_F0"][concStr] = []
                for pos in range(0, len(data[targets[tar]]["Cq"][concStr])):
                    data[targets[tar]]["Cq_F0"][concStr].append(np.log10(1.0 / np.power(data[targets[tar]]["cq_eff"], data[targets[tar]]["Cq"][concStr][pos])))
                    cq_check_y.append(np.log10(1.0 / np.power(data[targets[tar]]["cq_eff"], data[targets[tar]]["Cq"][concStr][pos])))
                data[targets[tar]]["cq_ss_per_dil"][concStr] = np.var(data[targets[tar]]["Cq_F0"][concStr], ddof=1) * (len(data[targets[tar]]["Cq_F0"][concStr]) - 1)
                data[targets[tar]]["cq_SSwithin"] += data[targets[tar]]["cq_ss_per_dil"][concStr]

            slope, intercept, r_val, p_val, std_err = scp.linregress(data[targets[tar]]["all_std_scaled_log"], cq_check_y)
            data[targets[tar]]["cq_slope_check"] = slope
            data[targets[tar]]["cq_ms_within_cq"] = data[targets[tar]]["cq_SSwithin"] / (data[targets[tar]]["var_n"] - data[targets[tar]]["var_groups"])
            data[targets[tar]]["cq_ms_within_linregpcr"] = data[targets[tar]]["SSwithin"] / (data[targets[tar]]["var_n"] - data[targets[tar]]["var_groups"])
            data[targets[tar]]["cq_ms_ratio"] = data[targets[tar]]["cq_ms_within_linregpcr"] / data[targets[tar]]["cq_ms_within_cq"]

        # difference calc
        dif_alpha = 0.95
        dif_detect_num = 0
        dif_detect_sum = 0.0
        data[targets[tar]]["dif_se_yfit"]  = {}
        data[targets[tar]]["dif_log_upper"]  = {}
        data[targets[tar]]["dif_log_lower"]  = {}
        data[targets[tar]]["dif_fold_up"]  = {}
        data[targets[tar]]["dif_fold_down"]  = {}
        dif_SSx = np.var(np.log10(data[targets[tar]]["all_std"]), ddof=1) * (data[targets[tar]]["var_n"] - 2.0)
        for concPos in range(0, len(data[targets[tar]]["std"]["val"])):
            concStr = data[targets[tar]]["std"]["str"][concPos]
            concFloat = data[targets[tar]]["std"]["val"][concPos]
            dif_x_mean = np.mean(np.log10(data[targets[tar]]["std"]["val"]))
            dif_n = len(data[targets[tar]]["scaled_log"][concStr])  # Any n of samples with same conc
            dif_t = -1.0 * scp.t.ppf((1.0 - dif_alpha) / 4.0, dif_n - 1)
            data[targets[tar]]["dif_se_yfit"][concStr] = data[targets[tar]]["stexy"] * np.power(1.0 / dif_n + np.power(np.log10(concFloat) - dif_x_mean, 2.0) / dif_SSx , 0.5)
            data[targets[tar]]["dif_log_upper"][concStr] = data[targets[tar]]["mean_log"][concStr] + dif_t * data[targets[tar]]["dif_se_yfit"][concStr]
            data[targets[tar]]["dif_log_lower"][concStr] = data[targets[tar]]["mean_log"][concStr] - dif_t * data[targets[tar]]["dif_se_yfit"][concStr]
            data[targets[tar]]["dif_fold_up"][concStr] = np.power(10.0, data[targets[tar]]["dif_log_upper"][concStr]) / np.power(10.0, data[targets[tar]]["mean_log"][concStr])
            dif_detect_sum += np.log(data[targets[tar]]["dif_fold_up"][concStr])
            dif_detect_num += 1
            data[targets[tar]]["dif_fold_down"][concStr] = np.power(10.0, data[targets[tar]]["mean_log"][concStr]) / np.power(10.0, data[targets[tar]]["dif_log_lower"][concStr])

        data[targets[tar]]["dif_detectable_diff"] = np.exp(dif_detect_sum / dif_detect_num)

    return data

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

        data = "<rdml version='1.4' xmlns:rdml='http://www.rdml.org' xmlns='http://www.rdml.org'>\n<dateMade>"
        data += datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
        data += "</dateMade>\n<dateUpdated>"
        data += datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
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

    def load_any_zip(self, filename):
        """Load an RDML file with decompression of first file. Uses loadXMLString().

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        if zipfile.is_zipfile(filename):
            self._rdmlFilename = filename
            zf = zipfile.ZipFile(filename, 'r')
            archiv_name = ""
            zip_list = zf.infolist()
            for curr_file in zip_list:
                if not curr_file.is_dir():
                    archiv_name = curr_file.filename
                    break
            if archiv_name != "":
                try:
                    data = zf.read(archiv_name).decode('utf-8')
                except KeyError:
                    raise RdmlError('No readable data in compressed RDML file found.')
                except UnicodeDecodeError:
                    raise RdmlError('No readable unicode data in compressed RDML file found.')
                else:
                    self.loadXMLString(data)
            zf.close()
            if archiv_name == "":
                raise RdmlError('No readable data in compressed RDML file found.')
        else:
            raise RdmlError('File format error, no compressed RDML file found.')

    def save(self, filename):
        """Save an RDML file with compression of rdml_data.xml.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to save to.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        elem = _get_or_create_subelement(self._node, "dateUpdated", self.xmlkeys())
        elem.text = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
        data = et.tostring(self._rdmlData, pretty_print=True)
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
            self._rdmlData = et.ElementTree(et.fromstring(data.encode('utf-8')))
            # Change to bytecode and defused?
        except et.XMLSyntaxError:
            raise RdmlError('XML load error, not a valid RDML or XML file.')
        self._node = self._rdmlData.getroot()
        if self._node.tag.replace("{http://www.rdml.org}", "") != 'rdml':
            raise RdmlError('Root element is not \'rdml\', not a valid RDML or XML file.')
        rdml_version = self._node.get('version')
        # Remainder: Update version in new() and validate()
        if rdml_version not in ['1.0', '1.1', '1.2', '1.3', '1.4']:
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
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        elif version == '1.3':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_3_REC.xsd'))
        elif version == '1.4':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_4_CR.xsd'))
        else:
            notes += 'RDML version:\tFalse\tUnknown schema version' + version + '\n'
            return notes
        notes += "RDML version:\tTrue\t" + version + "\n"

        xmlschema = et.XMLSchema(xmlschema_doc)
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
            except RdmlError:
                return False
        else:
            vd = self
        version = vd.version()
        rdmlws = os.path.dirname(os.path.abspath(__file__))
        if version == '1.0':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        elif version == '1.3':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_3_REC.xsd'))
        elif version == '1.4':
            xmlschema_doc = et.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_4_CR.xsd'))
        else:
            return False
        xmlschema = et.XMLSchema(xmlschema_doc)
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
                    "doubleStranded", "quantity", "calibratorSample", "cdnaSynthesisMethod",
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
                new_node = et.Element("dye", id=dkey)
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
                        old_letter = ord(re.sub(r"\d", "", old_id).upper()) - ord("A")
                        old_nr = int(re.sub(r"\D", "", old_id))
                        newId = old_nr + old_letter * int(columns)
                        node3.attrib['id'] = str(newId)
                if old_format == "3072-well plate; A1a1-D12h8":
                    exp3 = _get_all_children(node2, "react")
                    for node3 in exp3:
                        old_id = node3.get('id')
                        old_left = re.sub(r"\D\d+$", "", old_id)
                        old_left_letter = ord(re.sub(r"\d", "", old_left).upper()) - ord("A")
                        old_left_nr = int(re.sub(r"\D", "", old_left)) - 1
                        old_right = re.sub(r"^\D\d+", "", old_id)
                        old_right_letter = ord(re.sub(r"\d", "", old_right).upper()) - ord("A")
                        old_right_nr = int(re.sub(r"\D", "", old_right))
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
        hint2 = ""
        hint3 = ""
        hint4 = ""
        hint5 = ""
        hint6 = ""
        hint7 = ""
        hint8 = ""
        hint9 = ""
        hint10 = ""
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

                    exp4b = _get_all_children(node3, "data")
                    for node4 in exp4b:
                        exp5 = _get_all_children(node4, "ampEffMet")
                        for node5 in exp5:
                            hint2 = "Migration to v1.2 deleted react data \"ampEffMet\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "N0")
                        for node5 in exp5:
                            hint3 = "Migration to v1.2 deleted react data \"N0\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "ampEff")
                        for node5 in exp5:
                            hint4 = "Migration to v1.2 deleted react data \"ampEff\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "ampEffSE")
                        for node5 in exp5:
                            hint5 = "Migration to v1.2 deleted react data \"ampEffSE\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "corrF")
                        for node5 in exp5:
                            hint6 = "Migration to v1.2 deleted react data \"corrF\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "corrP")
                        for node5 in exp5:
                            hint7 = "Migration to v1.2 deleted react data \"corrP\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "corrCq")
                        for node5 in exp5:
                            hint8 = "Migration to v1.2 deleted react data \"corrCq\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "meltTemp")
                        for node5 in exp5:
                            hint9 = "Migration to v1.2 deleted react data \"meltTemp\" elements."
                            node4.remove(node5)
                        exp5 = _get_all_children(node4, "note")
                        for node5 in exp5:
                            hint10 = "Migration to v1.2 deleted react data \"note\" elements."
                            node4.remove(node5)
        if hint != "":
            ret.append(hint)
        if hint2 != "":
            ret.append(hint2)
        if hint3 != "":
            ret.append(hint3)
        if hint4 != "":
            ret.append(hint4)
        if hint5 != "":
            ret.append(hint5)
        if hint6 != "":
            ret.append(hint6)
        if hint7 != "":
            ret.append(hint7)
        if hint8 != "":
            ret.append(hint8)
        if hint9 != "":
            ret.append(hint9)
        if hint10 != "":
            ret.append(hint10)

        exp1 = _get_all_children(self._node, "sample")
        hint = ""
        hint2 = ""
        hint3 = ""
        hint4 = ""
        for node1 in exp1:
            exp2 = _get_all_children(node1, "type")
            if len(exp2) > 0:
                if "targetId" in exp2[0].attrib:
                    del exp2[0].attrib["targetId"]
                    hint = "Migration to v1.2 deleted sample type \"targetId\" attribute."
                for elCount in range(1, len(exp2)):
                    node1.remove(exp2[elCount])
                    hint2 = "Migration to v1.2 deleted sample \"type\" elements."
            exp2q = _get_all_children(node1, "quantity")
            if len(exp2q) > 0:
                if "targetId" in exp2q[0].attrib:
                    del exp2q[0].attrib["targetId"]
                    hint3 = "Migration to v1.2 deleted sample quantity \"targetId\" attribute."
                for elCount in range(1, len(exp2q)):
                    node1.remove(exp2q[elCount])
                    hint4 = "Migration to v1.2 deleted sample \"quantity\" elements."
        if hint != "":
            ret.append(hint)
        if hint2 != "":
            ret.append(hint2)
        if hint3 != "":
            ret.append(hint3)
        if hint4 != "":
            ret.append(hint4)

        exp1 = _get_all_children(self._node, "target")
        hint = ""
        for node1 in exp1:
            exp2 = _get_all_children(node1, "meltingTemperature")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted target \"meltingTemperature\" elements."
        if hint != "":
            ret.append(hint)

        exp1 = _get_all_children(self._node, "dye")
        hint = ""
        for node1 in exp1:
            exp2 = _get_all_children(node1, "dyeChemistry")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.2 deleted dye \"dyeChemistry\" elements."
        if hint != "":
            ret.append(hint)

        self._node.attrib['version'] = "1.2"
        return ret

    def migrate_version_1_3_to_1_4(self):
        """Migrates the rdml version from v1.3 to v1.4.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        copyThreshold = _get_copies_threshold()
        dataXMLelements = _getXMLDataType()
        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.3':
            raise RdmlError('RDML version for migration has to be v1.3.')

        self._node.attrib['version'] = "1.4"
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "react")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "data")
                    for node4 in exp4:
                        n0Val = _get_first_child_text(node4, "N0")
                        if n0Val == "":
                            continue
                        try:
                            n0Val = float(n0Val)
                        except ValueError:
                            continue
                        if not math.isfinite(n0Val):
                            continue
                        if n0Val <= 0.0:
                            continue
                        threshold = _get_first_child_text(node4, "quantFluor")
                        if threshold == "":
                            continue
                        try:
                            threshold = float(threshold)
                        except ValueError:
                            continue
                        if not math.isfinite(threshold):
                            continue
                        if threshold <= 0.0:
                            continue
                        Ncopy = copyThreshold * n0Val / threshold
                        _change_subelement(node4, "Ncopy", dataXMLelements, "{:.4f}".format(Ncopy), True, "string")
        return ret

    def migrate_version_1_4_to_1_3(self):
        """Migrates the rdml version from v1.4 to v1.3.

        Args:
            self: The class self parameter.

        Returns:
            A list of strings with the modifications made.
        """

        ret = []
        rdml_version = self._node.get('version')
        if rdml_version != '1.4':
            raise RdmlError('RDML version for migration has to be v1.4.')

        hint = ""
        hint2 = ""
        exp1 = _get_all_children(self._node, "experiment")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "run")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "react")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "vol")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react \"vol\" elements."
                        node3.remove(node4)

                    exp4b = _get_all_children(node3, "data")
                    for node4 in exp4b:
                        exp5 = _get_all_children(node4, "Ncopy")
                        for node5 in exp5:
                            hint2 = "Migration to v1.3 deleted react data \"Ncopy\" elements."
                            node4.remove(node5)
        if hint != "":
            ret.append(hint)
        if hint2 != "":
            ret.append(hint2)

        hint = ""
        exp1 = _get_all_children(self._node, "sample")
        for node1 in exp1:
            hint = ""
            exp2 = _get_all_children(node1, "doubleStranded")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.3 deleted sample \"doubleStranded\" element."
        if hint != "":
            ret.append(hint)

        hint = ""
        exp1 = _get_all_children(self._node, "target")
        for node1 in exp1:
            exp2 = _get_all_children(node1, "sequences")
            for node2 in exp2:
                exp3 = _get_all_children(node2, "forwardPrimer")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "oligoConc")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react data \"oligoConc\" elements."
                        node3.remove(node4)
                exp3 = _get_all_children(node2, "reversePrimer")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "oligoConc")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react data \"oligoConc\" elements."
                        node3.remove(node4)
                exp3 = _get_all_children(node2, "probe1")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "oligoConc")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react data \"oligoConc\" elements."
                        node3.remove(node4)
                exp3 = _get_all_children(node2, "probe2")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "oligoConc")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react data \"oligoConc\" elements."
                        node3.remove(node4)
                exp3 = _get_all_children(node2, "amplicon")
                for node3 in exp3:
                    exp4 = _get_all_children(node3, "oligoConc")
                    for node4 in exp4:
                        hint = "Migration to v1.3 deleted react data \"oligoConc\" elements."
                        node3.remove(node4)
        if hint != "":
            ret.append(hint)

        exp1 = _get_all_children(self._node, "dye")
        hint = ""
        hint2 = ""
        for node1 in exp1:
            exp2 = _get_all_children(node1, "dNTPs")
            for node2 in exp2:
                node1.remove(node2)
                hint = "Migration to v1.3 deleted dye \"dNTPs\" elements."
            exp2b = _get_all_children(node1, "dyeConc")
            for node2 in exp2b:
                node1.remove(node2)
                hint = "Migration to v1.3 deleted dye \"dyeConc\" elements."
        if hint != "":
            ret.append(hint)
        if hint2 != "":
            ret.append(hint2)

        self._node.attrib['version'] = "1.3"
        return ret

    def recreate_lost_ids(self):
        """Searches for lost ids and repairs them.

        Args:
            self: The class self parameter.

        Returns:
            A string with the modifications.
        """

        mess = ""
        rdml_version = self._node.get('version')

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
                forId = _get_first_child(subNode, "thermalCyclingConditions")
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
                        if lastNode.get('id') == None:
                            lastNode.attrib['id'] = "__no_sample_provided__"
                            foundIds["__no_sample_provided__"] = 0
                        else:
                            foundIds[lastNode.attrib['id']] = 0
        presentIds = []
        exp = _get_all_children(self._node, "sample")
        for node in exp:
            presentIds.append(node.attrib['id'])
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_sample(id=used_id, newposition=0)
                if rdml_version != '1.3':
                    elem = self.get_sample(byid=used_id)
                    elem.new_type("unkn", targetId=None)
                mess += "Recreated sample: " + used_id + "\n"
        # Find lost target
        foundIds = {}
        allExp = _get_all_children(self._node, "sample")
        for node in allExp:
            subNodes = _get_all_children(node, "type")
            for subNode in subNodes:
                if "targetId" in subNode.attrib:
                    foundIds[subNode.attrib['targetId']] = 0
            subNodes = _get_all_children(node, "quantity")
            for subNode in subNodes:
                if "targetId" in subNode.attrib:
                    foundIds[subNode.attrib['targetId']] = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    for dataNode in dataNodes:
                        noTarget = True
                        lastNodes = _get_all_children(dataNode, "tar")
                        for lastNode in lastNodes:
                            if lastNode.get('id') == None:
                                lastNode.attrib['id'] = "__no_target_provided__"
                                foundIds["__no_target_provided__"] = 0
                            else:
                                foundIds[lastNode.attrib['id']] = 0
                            noTarget = False
                        if noTarget:
                            new_node = et.Element("tar")
                            new_node.attrib['id'] = "__no_target_provided__"
                            dataNode.insert(0, new_node)
                            foundIds["__no_target_provided__"] = 0
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
        recreateDye = False
        exp = _get_all_children(self._node, "target")
        for node in exp:
            presentIds.append(node.attrib['id'])
        rdTars = self.targets()
        for rdTar in rdTars:
            tarType = _get_first_child_text(rdTar._node, "type")
            if tarType == "":
                rdTar["type"] = "toi"
        for used_id in foundIds:
            if used_id not in presentIds:
                self.new_target(id=used_id, type="toi", newposition=0)
                elem = self.get_target(byid=used_id)
                elem["dyeId"] = "recreated dye"
                mess += "Recreated target: " + used_id + "\n"
                recreateDye = True
        presentIds = []
        if recreateDye:
            exp = _get_all_children(self._node, "dye")
            for node in exp:
                presentIds.append(node.attrib['id'])
            if "recreated dye" not in presentIds:
                self.new_dye(id="recreated dye", newposition=0)
                mess += "Recreated new dye: recreated dye\n"
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
        count = 0
        allExp = _get_all_children(self._node, "experiment")
        for node in allExp:
            subNodes = _get_all_children(node, "run")
            for subNode in subNodes:
                foundIds = {}
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

    def fixTempsMeltcurve(self):
        """Searches in experiment-run-react for different temps in one read and keeps the mean.

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
                # Check we have always the same number of columns
                colCount = -10
                rowCount = 0
                reactNodes = _get_all_children(subNode, "react")
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    rowCount += 1
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "mdp")
                        if colCount < 0:
                            colCount = len(lastNodes)
                        else:
                            if colCount != len(lastNodes):
                                mess = "The wells have different number of temperatures in the "
                                mess += "melting curve, fixing is not possible!\n"
                                return mess
                # Read the temps in a matrix
                allTemps = np.zeros((rowCount, colCount), dtype=np.float64)
                currRow = -1
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    currRow += 1
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "mdp")
                        currCol = -1
                        for mdp in lastNodes:
                            currCol += 1
                            allTemps[currRow, currCol] = float(_get_first_child_text(mdp, "tmp"))
                meanTemp = np.mean(allTemps, axis=0)
                # Rewrite the temperatures
                currRow = -1
                for reactNode in reactNodes:
                    dataNodes = _get_all_children(reactNode, "data")
                    currRow += 1
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "mdp")
                        currCol = -1
                        for mdp in lastNodes:
                            currCol += 1
                            editNode = _get_first_child(mdp, "tmp")
                            editNode.text = str("%.2f" % float(meanTemp[currCol]))
                            count += 1
        if count > 0:
            mess = "The temperatures in the melting curve were " + str(count) + " times rewritten!\n"

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

        new_node = et.Element("id")
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

    def import_experimenter(self, experimenter):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            experimenter: The target to import as Experimenter element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), 999999)
        currId = experimenter["id"]
        allChildren = _get_all_children(self._node, "experimenter")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, experimenter._node)

    def import_all_experimenters(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the experimenters from
            addMode: "all" (default) , "only-new" or "update"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        for exp in self.experimenters():
            known[exp["id"]] = 1

        for addEle in add_rd.experimenters():
            if addMode == "only-new":
                if addEle["id"] not in known:
                    self.import_experimenter(addEle)
            elif addMode == "update":
                if addEle["id"] in known:
                    self.import_experimenter(addEle)
            else:
                self.import_experimenter(addEle)

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

    def import_documentation(self, documentation):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            documentation: The target to import as Experimenter element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999)
        currId = documentation["id"]
        allChildren = _get_all_children(self._node, "documentation")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, documentation._node)

    def import_all_documentations(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the documentations from
            addMode: "all" (default) , "only-new" or "update"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        for exp in self.documentations():
            known[exp["id"]] = 1

        for addEle in add_rd.documentations():
            if addMode == "only-new":
                if addEle["id"] not in known:
                    self.import_documentation(addEle)
            elif addMode == "update":
                if addEle["id"] in known:
                    self.import_documentation(addEle)
            else:
                self.import_documentation(addEle)

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

    def import_dye(self, dye):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            dye: The target to import as Experimenter element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "dye", self.xmlkeys(), 999999)
        currId = dye["id"]
        allChildren = _get_all_children(self._node, "dye")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "dye", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, dye._node)

    def import_all_dyes(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the dyes from
            addMode: "all" (default) , "only-new" or "update"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        for exp in self.dyes():
            known[exp["id"]] = 1

        for addEle in add_rd.dyes():
            if addMode == "only-new":
                if addEle["id"] not in known:
                    self.import_dye(addEle)
            elif addMode == "update":
                if addEle["id"] in known:
                    self.import_dye(addEle)
            else:
                self.import_dye(addEle)

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

    def new_sample(self, id, newposition=None):
        """Creates a new sample element.

        Args:
            self: The class self parameter.
            id: Sample unique id (required)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = _create_new_element(self._node, "sample", id)
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

    def import_sample(self, add_rd, sample, addMode):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the samples from
            sample: The sample to import as Experimenter element
            addMode: "only-new" (default) , "all-dep" or "no-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if addMode != "no-dep":
            presDocs = []
            addDocs = []
            presDyes = []
            addDyes = []
            presTargets = []
            addTargets = []
            presThermal = []
            addThermal = []
            presExps = []
            addExps = []

            exp = _get_all_children(self._node, "documentation")
            for node in exp:
                presDocs.append(node.attrib['id'])
            exp = _get_all_children(self._node, "dye")
            for node in exp:
                presDyes.append(node.attrib['id'])
            exp = _get_all_children(self._node, "target")
            for node in exp:
                presTargets.append(node.attrib['id'])
            exp = _get_all_children(self._node, "thermalCyclingConditions")
            for node in exp:
                presThermal.append(node.attrib['id'])
            exp = _get_all_children(self._node, "experimenter")
            for node in exp:
                presExps.append(node.attrib['id'])

            subNodes = _get_all_children(sample._node, "documentation")
            for subNode in subNodes:
                foundId = subNode.attrib['id']
                if foundId not in addDocs:
                    if foundId not in presDocs:
                        addDocs.append(foundId)
                    if foundId not in addDocs:
                        if addMode == "all-dep":
                            addDocs.append(foundId)

            subNodes = _get_all_children(sample._node, "type")
            for subNode in subNodes:
                if 'targetId' in subNode.attrib:
                    foundId = subNode.attrib['targetId']
                    if foundId not in addTargets:
                        if foundId not in presTargets:
                            addTargets.append(foundId)
                        if foundId not in addTargets:
                            if addMode == "all-dep":
                                addTargets.append(foundId)

            subNodes = _get_all_children(sample._node, "quantity")
            for subNode in subNodes:
                if 'targetId' in subNode.attrib:
                    foundId = subNode.attrib['targetId']
                    if foundId not in addTargets:
                        if foundId not in presTargets:
                            addTargets.append(foundId)
                        if foundId not in addTargets:
                            if addMode == "all-dep":
                                addTargets.append(foundId)

            subNode = _get_first_child(sample._node, "cdnaSynthesisMethod")
            if subNode is not None:
                forId = _get_first_child(subNode, "thermalCyclingConditions")
                if forId is not None:
                    foundId = forId.attrib['id']
                    if foundId not in addThermal:
                        if foundId not in presThermal:
                            addThermal.append(foundId)
                        if foundId not in addThermal:
                            if addMode == "all-dep":
                                addThermal.append(foundId)

            # Get the dependent of the dependent
            exp = _get_all_children(add_rd._node, "target")
            for node in exp:
                if node.attrib['id'] in addTargets:
                    subNodes = _get_all_children(node, "experimenter")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addExps:
                            if foundId not in presExps:
                                addExps.append(foundId)
                            if foundId not in addExps:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addExps.append(foundId)
                    subNode = _get_first_child(node, "dyeId")
                    if subNode is not None:
                        foundId = subNode.attrib['id']
                        if foundId not in addDyes:
                            if foundId not in presDyes:
                                addDyes.append(foundId)
                            if foundId not in addDyes:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addDyes.append(foundId)
            exp = _get_all_children(add_rd._node, "thermalCyclingConditions")
            for node in exp:
                if node.attrib['id'] in addThermal:
                    subNodes = _get_all_children(node, "experimenter")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addExps:
                            if foundId not in presExps:
                                addExps.append(foundId)
                            if foundId not in addExps:
                                if addMode == "all-dep":
                                    addExps.append(foundId)
                    subNodes = _get_all_children(node, "documentation")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addDocs:
                            if foundId not in presDocs:
                                addDocs.append(foundId)
                            if foundId not in addDocs:
                                if addMode == "all-dep":
                                    addDocs.append(foundId)

            for docId in addDocs:
                self.import_documentation(add_rd.get_documentation(byid=docId))
            for docId in addDyes:
                self.import_dye(add_rd.get_dye(byid=docId))
            for docId in addTargets:
                self.import_target(add_rd, add_rd.get_target(byid=docId), "no-dep")
            for docId in addThermal:
                self.import_therm_cyc_cons(add_rd, add_rd.get_therm_cyc_cons(byid=docId), "no-dep")
            for docId in addExps:
                self.import_experimenter(add_rd.get_experimenter(byid=docId))

        pos = _get_tag_pos(self._node, "sample", self.xmlkeys(), 999999)
        currId = sample["id"]
        allChildren = _get_all_children(self._node, "sample")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "sample", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, sample._node)

    def import_all_samples(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the dyes from
            addMode: "all" (default) , "only-new", "update", "all-incl-dep",
                     "only-new-incl-dep" or "update-incl-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        presDocs = []
        addDocs = []
        presDyes = []
        addDyes = []
        presTargets = []
        addTargets = []
        presThermal = []
        addThermal = []
        presExps = []
        addExps = []

        for exp in self.samples():
            known[exp["id"]] = 1

        exp = _get_all_children(self._node, "documentation")
        for node in exp:
            presDocs.append(node.attrib['id'])
        exp = _get_all_children(self._node, "dye")
        for node in exp:
            presDyes.append(node.attrib['id'])
        exp = _get_all_children(self._node, "target")
        for node in exp:
            presTargets.append(node.attrib['id'])
        exp = _get_all_children(self._node, "thermalCyclingConditions")
        for node in exp:
            presThermal.append(node.attrib['id'])
        exp = _get_all_children(self._node, "experimenter")
        for node in exp:
            presExps.append(node.attrib['id'])

        allNew = add_rd.samples()
        for addEle in allNew:
            readChild = False
            if addMode in ["only-new", "only-new-incl-dep"]:
                if addEle["id"] not in known:
                    self.import_sample(add_rd, addEle, "no-dep")
                    readChild = True
            elif addMode in ["update", "update-incl-dep"]:
                if addEle["id"] in known:
                    self.import_sample(add_rd, addEle, "no-dep")
                    readChild = True
            else:
                self.import_sample(add_rd, addEle, "no-dep")
                readChild = True

            if readChild:
                subNodes = _get_all_children(addEle._node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)

                subNodes = _get_all_children(addEle._node, "type")
                for subNode in subNodes:
                    if 'targetId' in subNode.attrib:
                        foundId = subNode.attrib['targetId']
                        if foundId not in addTargets:
                            if foundId not in presTargets:
                                addTargets.append(foundId)
                            if foundId not in addTargets:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addTargets.append(foundId)

                subNodes = _get_all_children(addEle._node, "quantity")
                for subNode in subNodes:
                    if 'targetId' in subNode.attrib:
                        foundId = subNode.attrib['targetId']
                        if foundId not in addTargets:
                            if foundId not in presTargets:
                                addTargets.append(foundId)
                            if foundId not in addTargets:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addTargets.append(foundId)

                subNode = _get_first_child(addEle._node, "cdnaSynthesisMethod")
                if subNode is not None:
                    forId = _get_first_child(subNode, "thermalCyclingConditions")
                    if forId is not None:
                        foundId = forId.attrib['id']
                        if foundId not in addThermal:
                            if foundId not in presThermal:
                                addThermal.append(foundId)
                            if foundId not in addThermal:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addThermal.append(foundId)

        # Get the dependent of the dependent
        exp = _get_all_children(add_rd._node, "target")
        for node in exp:
            if node.attrib['id'] in addTargets:
                subNodes = _get_all_children(node, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addExps.append(foundId)
                subNode = _get_first_child(node, "dyeId")
                if subNode is not None:
                    foundId = subNode.attrib['id']
                    if foundId not in addDyes:
                        if foundId not in presDyes:
                            addDyes.append(foundId)
                        if foundId not in addDyes:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDyes.append(foundId)
        exp = _get_all_children(add_rd._node, "thermalCyclingConditions")
        for node in exp:
            if node.attrib['id'] in addThermal:
                subNodes = _get_all_children(node, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addExps.append(foundId)
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)

        for docId in addDocs:
            self.import_documentation(add_rd.get_documentation(byid=docId))
        for docId in addDyes:
            self.import_dye(add_rd.get_dye(byid=docId))
        for docId in addTargets:
            self.import_target(add_rd, add_rd.get_target(byid=docId), "no-dep")
        for docId in addThermal:
            self.import_therm_cyc_cons(add_rd, add_rd.get_therm_cyc_cons(byid=docId), "no-dep")
        for docId in addExps:
            self.import_experimenter(add_rd.get_experimenter(byid=docId))

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

    def export_annotations(self):
        """Returns a table the samples with the annotations in a string.

        Args:
            self: The class self parameter.

        Returns:
            A table with \t seperators.
        """

        ret = "Sample"
        lookupAnno = {}

        # Collect the unique annotation keys
        exp = _get_all_children(self._node, "sample")
        for node in exp:
            annos = _get_all_children(node, "annotation")
            for node in annos:
                prop = _get_first_child_text(node, "property")
                if prop != "":
                    lookupAnno[prop] = 1

        # Sort the keys and correct lookup
        annoSort = list(lookupAnno.keys())
        annoSort.sort()
        anno_pos = 1
        for anno in annoSort:
            ret += '\t' + anno
            lookupAnno[anno] = anno_pos
            anno_pos += 1
        ret += '\n'

        # Fill the table grid
        resTab = []
        sam_pos = 0
        for node in exp:
            resTab.append([node.attrib['id']])
            for anno in annoSort:
                resTab[sam_pos].append("")
            annos = _get_all_children(node, "annotation")
            for node in annos:
                prop = _get_first_child_text(node, "property")
                value = _get_first_child_text(node, "value")
                if prop != "":
                    resTab[sam_pos][lookupAnno[prop]] = value
            sam_pos += 1

        # Print Table out
        for row in resTab:
            for colTxt in row:
                ret += colTxt + '\t' 
            ret = re.sub(r"\t$", "\n", ret)

        return ret

    def import_annotations(self, csvData):
        """Returns a list of the annotations in the xml file.

        Args:
            self: The class self parameter.
            csvData: A tsv file with the annotation data

        Returns:
            Nothing.
        """

        ver = self._node.get('version')
        if ver == "1.1":
            return "Error: Not available in RDML v1.1."

        err = ""
        presSamples = []
        exp = _get_all_children(self._node, "sample")
        for node in exp:
            presSamples.append(node.attrib['id'])
        with open(csvData, newline='') as tfile:  # add encoding='utf-8' ?
            try:
                annoTab = list(csv.reader(tfile, delimiter='\t'))
            except UnicodeDecodeError:
                raise RdmlError('The annotation must have UTF-8 encoding.')
            if len(annoTab) < 2:
                raise RdmlError('The annotation file must have at least two rows.')
            if len(annoTab[0]) < 2:
                raise RdmlError('The annotation file must have at least two columns.')
            for row in range(1, len(annoTab)):
                samId = annoTab[row][0]
                if samId != "":
                    if samId not in presSamples:
                        if err == "":
                            err += samId
                        else:
                            err += ", " + samId
                        continue
                    for col in range(1, len(annoTab[0])):
                        annoProp = annoTab[0][col]
                        annoVal = annoTab[row][col]
                        if annoProp != "":
                            el = self.get_sample(byid=samId)
                            el.edit_annotation_value(property=annoProp, value=annoVal)
        if err != "":
            err = "Error: The uploaded sample ids are not found in the RDML file: " + err
        return err

    def rename_annotation_property(self, oldProperty, newProperty):
        """Returns a list of the annotations in the xml file.

        Args:
            self: The class self parameter.
            oldProperty: The old property value to be replaced
            newProperty: The new property value to use

        Returns:
            A list of dics with property and value strings.
        """

        ver = self._node.get('version')
        if ver == "1.1":
            return ""

        exp = _get_all_children(self._node, "sample")
        for node in exp:
            annos = _get_all_children(node, "annotation")
            for node2 in annos:
                prop = _get_first_child(node2, "property")
                if prop.text == oldProperty:
                    prop.text = newProperty
        return

    def rename_annotation_value(self, property, oldValue, newValue):
        """Returns a list of the annotations in the xml file.

        Args:
            self: The class self parameter.
            property: The property to work in
            oldValue: The old value to be replaced
            newProperty: The new value to use

        Returns:
            A list of dics with property and value strings.
        """

        ver = self._node.get('version')
        if ver == "1.1":
            return ""

        exp = _get_all_children(self._node, "sample")
        for node in exp:
            annos = _get_all_children(node, "annotation")
            for node2 in annos:
                prop = _get_first_child(node2, "property")
                value = _get_first_child(node2, "value")
                if prop.text == property:
                    if value.text == oldValue:
                        value.text = newValue
        return

    def combine_annotations(self, combinedPoperty, leftProperty, connectProperty, rightProperty):
        """Returns a list of the annotations in the xml file.

        Args:
            self: The class self parameter.
            combinedPoperty: The new property name
            leftProperty: The property to take the left part of the new value
            connectProperty: A string to place between the new values
            rightProperty: The property to take the right part of the new value

        Returns:
            A list of dics with property and value strings.
        """

        ver = self._node.get('version')
        if ver == "1.1":
            return ""

        exp = _get_all_children(self._node, "sample")
        for node in exp:
            leftVal = ""
            rightVal = ""
            annos = _get_all_children(node, "annotation")
            for node2 in annos:
                prop = _get_first_child_text(node2, "property")
                if prop == leftProperty:
                    leftVal = _get_first_child_text(node2, "value")
                if prop == rightProperty:
                    rightVal = _get_first_child_text(node2, "value")
            if leftVal != "":
                if rightVal != "":
                    el = Sample(node)
                    newVal = leftVal + connectProperty + rightVal
                    el.edit_annotation_value(property=combinedPoperty, value=newVal)
        return

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

    def import_target(self, add_rd, target, addMode):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the targets from
            target: The target to import as Experimenter element
            addMode: "only-new" (default) , "all-dep" or "no-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if addMode != "no-dep":
            presDocs = []
            addDocs = []
            presDyes = []
            addDyes = []

            exp = _get_all_children(self._node, "documentation")
            for node in exp:
                presDocs.append(node.attrib['id'])
            exp = _get_all_children(self._node, "dye")
            for node in exp:
                presDyes.append(node.attrib['id'])

            subNodes = _get_all_children(target._node, "documentation")
            for subNode in subNodes:
                foundId = subNode.attrib['id']
                if foundId not in addDocs:
                    if foundId not in presDocs:
                        addDocs.append(foundId)
                    if foundId not in addDocs:
                        if addMode == "all-dep":
                            addDocs.append(foundId)
            subNode = _get_first_child(target._node, "dyeId")
            if subNode is not None:
                foundId = subNode.attrib['id']
                if foundId not in addDyes:
                    if foundId not in presDyes:
                        addDyes.append(foundId)
                    if foundId not in addDyes:
                        if addMode == "all-dep":
                            addDyes.append(foundId)

            for docId in addDocs:
                self.import_documentation(add_rd.get_documentation(byid=docId))
            for docId in addDyes:
                self.import_dye(add_rd.get_dye(byid=docId))

        pos = _get_tag_pos(self._node, "target", self.xmlkeys(), 999999)
        currId = target["id"]
        allChildren = _get_all_children(self._node, "target")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "target", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, target._node)

    def import_all_targets(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the dyes from
            addMode: "all" (default) , "only-new", "update", "all-incl-dep",
                     "only-new-incl-dep" or "update-incl-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        presDocs = []
        addDocs = []
        presDyes = []
        addDyes = []
        for exp in self.targets():
            known[exp["id"]] = 1

        exp = _get_all_children(self._node, "documentation")
        for node in exp:
            presDocs.append(node.attrib['id'])

        exp = _get_all_children(self._node, "dye")
        for node in exp:
            presDyes.append(node.attrib['id'])

        allNew = add_rd.targets()
        for addEle in allNew:
            readChild = False
            if addMode in ["only-new", "only-new-incl-dep"]:
                if addEle["id"] not in known:
                    self.import_target(add_rd, addEle, "no-dep")
                    readChild = True
            elif addMode in ["update", "update-incl-dep"]:
                if addEle["id"] in known:
                    self.import_target(add_rd, addEle, "no-dep")
                    readChild = True
            else:
                self.import_target(add_rd, addEle, "no-dep")
                readChild = True

            if readChild:
                subNodes = _get_all_children(addEle._node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)
                subNode = _get_first_child(addEle._node, "dyeId")
                if subNode is not None:
                    foundId = subNode.attrib['id']
                    if foundId not in addDyes:
                        if foundId not in presDyes:
                            addDyes.append(foundId)
                        if foundId not in addDyes:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDyes.append(foundId)

        for docId in addDocs:
            self.import_documentation(add_rd.get_documentation(byid=docId))

        for docId in addDyes:
            self.import_dye(add_rd.get_dye(byid=docId))

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
        step = et.SubElement(new_node, "step")
        et.SubElement(step, "nr").text = "1"
        et.SubElement(step, "lidOpen")
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

    def import_therm_cyc_cons(self, add_rd, thermCycCon, addMode):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the therm_cyc_cons from
            thermCycCon: The thermCycCon to import as Experimenter element
            addMode: "only-new" (default) , "all-dep" or "no-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if addMode != "no-dep":
            presDocs = []
            addDocs = []
            presExps = []
            addExps = []

            exp = _get_all_children(self._node, "documentation")
            for node in exp:
                presDocs.append(node.attrib['id'])
            exp = _get_all_children(self._node, "experimenter")
            for node in exp:
                presExps.append(node.attrib['id'])

            subNodes = _get_all_children(thermCycCon._node, "documentation")
            for subNode in subNodes:
                foundId = subNode.attrib['id']
                if foundId not in addDocs:
                    if foundId not in presDocs:
                        addDocs.append(foundId)
                    if foundId not in addDocs:
                        if addMode == "all-dep":
                            addDocs.append(foundId)
            subNodes = _get_all_children(thermCycCon._node, "experimenter")
            for subNode in subNodes:
                foundId = subNode.attrib['id']
                if foundId not in addExps:
                    if foundId not in presExps:
                        addExps.append(foundId)
                    if foundId not in addExps:
                        if addMode == "all-dep":
                            addExps.append(foundId)

            for docId in addDocs:
                self.import_documentation(add_rd.get_documentation(byid=docId))
            for docId in addExps:
                self.import_experimenter(add_rd.get_experimenter(byid=docId))

        pos = _get_tag_pos(self._node, "thermalCyclingConditions", self.xmlkeys(), 999999)
        currId = thermCycCon["id"]
        allChildren = _get_all_children(self._node, "thermalCyclingConditions")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "thermalCyclingConditions", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, thermCycCon._node)

    def import_all_therm_cyc_cons(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the dyes from
            addMode: "all" (default) , "only-new", "update", "all-incl-dep",
                     "only-new-incl-dep" or "update-incl-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        presDocs = []
        addDocs = []
        presExps = []
        addExps = []
        for exp in self.therm_cyc_cons():
            known[exp["id"]] = 1

        exp = _get_all_children(self._node, "documentation")
        for node in exp:
            presDocs.append(node.attrib['id'])
        exp = _get_all_children(self._node, "experimenter")
        for node in exp:
            presExps.append(node.attrib['id'])

        allNew = add_rd.therm_cyc_cons()
        for addEle in allNew:
            readChild = False
            if addMode in ["only-new", "only-new-incl-dep"]:
                if addEle["id"] not in known:
                    self.import_therm_cyc_cons(add_rd, addEle, "no-dep")
                    readChild = True
            elif addMode in ["update", "update-incl-dep"]:
                if addEle["id"] in known:
                    self.import_therm_cyc_cons(add_rd, addEle, "no-dep")
                    readChild = True
            else:
                self.import_therm_cyc_cons(add_rd, addEle, "no-dep")
                readChild = True

            if readChild:
                subNodes = _get_all_children(addEle._node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)
                subNodes = _get_all_children(addEle._node, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addExps.append(foundId)

        for docId in addDocs:
            self.import_documentation(add_rd.get_documentation(byid=docId))
        for docId in addExps:
            self.import_experimenter(add_rd.get_experimenter(byid=docId))

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

    def move_experiment_run(self, id, oldid, runid):
        """Moves a run to a different experiment.

        Args:
            self: The class self parameter.
            id: The experiments unique id to insert to
            oldid: The experiments unique id with the current run
            runid: The runs unique id

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        newExp = self.get_experiment(byid=id)
        oldExp = self.get_experiment(byid=oldid)
        theRun = oldExp.get_run(byid=runid)
        place = _get_tag_pos(newExp._node, "run", newExp.xmlkeys(), 999999999)
        newExp._node.insert(place, theRun._node)

    def import_experiment(self, add_rd, experiment, addMode):
        """Imports the element to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the experiments from
            experiment: The experiment to import as Experimenter element
            addMode: "only-new" (default) , "all-dep" or "no-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if addMode != "no-dep":
            presExps = []
            addExps = []
            presDocs = []
            addDocs = []
            presDyes = []
            addDyes = []
            presSamples = []
            addSamples = []
            presTargets = []
            addTargets = []
            presThermal = []
            addThermal = []

            exp = _get_all_children(self._node, "experimenter")
            for node in exp:
                presExps.append(node.attrib['id'])
            exp = _get_all_children(self._node, "documentation")
            for node in exp:
                presDocs.append(node.attrib['id'])
            exp = _get_all_children(self._node, "dye")
            for node in exp:
                presDyes.append(node.attrib['id'])
            exp = _get_all_children(self._node, "sample")
            for node in exp:
                presSamples.append(node.attrib['id'])
            exp = _get_all_children(self._node, "target")
            for node in exp:
                presTargets.append(node.attrib['id'])
            exp = _get_all_children(self._node, "thermalCyclingConditions")
            for node in exp:
                presThermal.append(node.attrib['id'])

            subNodes = _get_all_children(experiment._node, "documentation")
            for subNode in subNodes:
                foundId = subNode.attrib['id']
                if foundId not in addDocs:
                    if foundId not in presDocs:
                        addDocs.append(foundId)
                    if foundId not in addDocs:
                        if addMode == "all-dep":
                            addDocs.append(foundId)
            runNodes = _get_all_children(experiment._node, "run")
            for runNode in runNodes:
                subNodes = _get_all_children(runNode, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode == "all-dep":
                                addDocs.append(foundId)
                subNodes = _get_all_children(runNode, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode == "all-dep":
                                addExps.append(foundId)
                forId = _get_first_child(runNode, "thermalCyclingConditions")
                if forId is not None:
                    foundId = forId.attrib['id']
                    if foundId not in addThermal:
                        if foundId not in presThermal:
                            addThermal.append(foundId)
                        if foundId not in addThermal:
                            if addMode == "all-dep":
                                addThermal.append(foundId)
                reactNodes = _get_all_children(runNode, "react")
                for reactNode in reactNodes:
                    sampleNode = _get_first_child(reactNode, "sample")
                    if sampleNode is not None:
                        foundId = sampleNode.attrib['id']
                        if foundId not in addSamples:
                            if foundId not in presSamples:
                                addSamples.append(foundId)
                            if foundId not in addSamples:
                                if addMode == "all-dep":
                                    addSamples.append(foundId)
                    dataNodes = _get_all_children(reactNode, "data")
                    for dataNode in dataNodes:
                        lastNodes = _get_all_children(dataNode, "tar")
                        for lastNode in lastNodes:
                            foundId = lastNode.attrib['id']
                            if foundId not in addTargets:
                                if foundId not in presTargets:
                                    addTargets.append(foundId)
                                if foundId not in addTargets:
                                    if addMode == "all-dep":
                                        addTargets.append(foundId)
                    partNodes = _get_all_children(reactNode, "partitions")
                    for partNode in partNodes:
                        dataNodes = _get_all_children(partNode, "data")
                        for dataNode in dataNodes:
                            lastNodes = _get_all_children(dataNode, "tar")
                            for lastNode in lastNodes:
                                foundId = lastNode.attrib['id']
                                if foundId not in addTargets:
                                    if foundId not in presTargets:
                                        addTargets.append(foundId)
                                    if foundId not in addTargets:
                                        if addMode == "all-dep":
                                            addTargets.append(foundId)

            # Get the dependent of the dependent
            exp = _get_all_children(add_rd._node, "sample")
            for node in exp:
                if node.attrib['id'] in addSamples:
                    subNodes = _get_all_children(node, "documentation")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addDocs:
                            if foundId not in presDocs:
                                addDocs.append(foundId)
                            if foundId not in addDocs:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addDocs.append(foundId)
                    subNodes = _get_all_children(node, "type")
                    for subNode in subNodes:
                        if 'targetId' in subNode.attrib:
                            foundId = subNode.attrib['targetId']
                            if foundId not in addTargets:
                                if foundId not in presTargets:
                                    addTargets.append(foundId)
                                if foundId not in addTargets:
                                    if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                        addTargets.append(foundId)
                    subNodes = _get_all_children(node, "quantity")
                    for subNode in subNodes:
                        if 'targetId' in subNode.attrib:
                            foundId = subNode.attrib['targetId']
                            if foundId not in addTargets:
                                if foundId not in presTargets:
                                    addTargets.append(foundId)
                                if foundId not in addTargets:
                                    if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                        addTargets.append(foundId)
                    subNode = _get_first_child(node, "cdnaSynthesisMethod")
                    if subNode is not None:
                        forId = _get_first_child(subNode, "thermalCyclingConditions")
                        if forId is not None:
                            foundId = forId.attrib['id']
                            if foundId not in addThermal:
                                if foundId not in presThermal:
                                    addThermal.append(foundId)
                                if foundId not in addThermal:
                                    if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                        addThermal.append(foundId)
            exp = _get_all_children(add_rd._node, "target")
            for node in exp:
                if node.attrib['id'] in addTargets:
                    subNodes = _get_all_children(node, "experimenter")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addExps:
                            if foundId not in presExps:
                                addExps.append(foundId)
                            if foundId not in addExps:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addExps.append(foundId)
                    subNode = _get_first_child(node, "dyeId")
                    if subNode is not None:
                        foundId = subNode.attrib['id']
                        if foundId not in addDyes:
                            if foundId not in presDyes:
                                addDyes.append(foundId)
                            if foundId not in addDyes:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addDyes.append(foundId)
            exp = _get_all_children(add_rd._node, "thermalCyclingConditions")
            for node in exp:
                if node.attrib['id'] in addThermal:
                    subNodes = _get_all_children(node, "experimenter")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addExps:
                            if foundId not in presExps:
                                addExps.append(foundId)
                            if foundId not in addExps:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addExps.append(foundId)
                    subNodes = _get_all_children(node, "documentation")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addDocs:
                            if foundId not in presDocs:
                                addDocs.append(foundId)
                            if foundId not in addDocs:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addDocs.append(foundId)

            for docId in addExps:
                self.import_experimenter(add_rd.get_experimenter(byid=docId))
            for docId in addDocs:
                self.import_documentation(add_rd.get_documentation(byid=docId))
            for docId in addDyes:
                self.import_dye(add_rd.get_dye(byid=docId))
            for docId in addSamples:
                self.import_sample(add_rd, add_rd.get_sample(byid=docId), "no-dep")
            for docId in addTargets:
                self.import_target(add_rd, add_rd.get_target(byid=docId), "no-dep")
            for docId in addThermal:
                self.import_therm_cyc_cons(add_rd, add_rd.get_therm_cyc_cons(byid=docId), "no-dep")

        pos = _get_tag_pos(self._node, "experiment", self.xmlkeys(), 999999)
        currId = experiment["id"]
        allChildren = _get_all_children(self._node, "experiment")
        count = -1
        for node in allChildren:
            count += 1
            if node.get('id') == currId:
                pos = _get_tag_pos(self._node, "experiment", self.xmlkeys(), count)
                self._node.remove(node)
        self._node.insert(pos, experiment._node)

    def import_all_experiments(self, add_rd, addMode):
        """Imports all elements to the end of the current list.

        Args:
            self: The class self parameter.
            add_rd: The rdml to import the dyes from
            addMode: "all" (default) , "only-new", "update", "all-incl-dep",
                     "only-new-incl-dep" or "update-incl-dep"

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        known = {}
        presExps = []
        addExps = []
        presDocs = []
        addDocs = []
        presDyes = []
        addDyes = []
        presSamples = []
        addSamples = []
        presTargets = []
        addTargets = []
        presThermal = []
        addThermal = []

        for exp in self.experiments():
            known[exp["id"]] = 1

        exp = _get_all_children(self._node, "experimenter")
        for node in exp:
            presExps.append(node.attrib['id'])
        exp = _get_all_children(self._node, "documentation")
        for node in exp:
            presDocs.append(node.attrib['id'])
        exp = _get_all_children(self._node, "dye")
        for node in exp:
            presDyes.append(node.attrib['id'])
        exp = _get_all_children(self._node, "sample")
        for node in exp:
            presSamples.append(node.attrib['id'])
        exp = _get_all_children(self._node, "target")
        for node in exp:
            presTargets.append(node.attrib['id'])
        exp = _get_all_children(self._node, "thermalCyclingConditions")
        for node in exp:
            presThermal.append(node.attrib['id'])

        allNew = add_rd.experiments()
        for addEle in allNew:
            readChild = False
            if addMode in ["only-new", "only-new-incl-dep"]:
                if addEle["id"] not in known:
                    self.import_experiment(add_rd, addEle, "no-dep")
                    readChild = True
            elif addMode in ["update", "update-incl-dep"]:
                if addEle["id"] in known:
                    self.import_experiment(add_rd, addEle, "no-dep")
                    readChild = True
            else:
                self.import_experiment(add_rd, addEle, "no-dep")
                readChild = True

            if readChild:
                subNodes = _get_all_children(addEle._node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)
                runNodes = _get_all_children(addEle._node, "run")
                for runNode in runNodes:
                    subNodes = _get_all_children(runNode, "documentation")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addDocs:
                            if foundId not in presDocs:
                                addDocs.append(foundId)
                            if foundId not in addDocs:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addDocs.append(foundId)
                    subNodes = _get_all_children(runNode, "experimenter")
                    for subNode in subNodes:
                        foundId = subNode.attrib['id']
                        if foundId not in addExps:
                            if foundId not in presExps:
                                addExps.append(foundId)
                            if foundId not in addExps:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addExps.append(foundId)
                    forId = _get_first_child(runNode, "thermalCyclingConditions")
                    if forId is not None:
                        foundId = forId.attrib['id']
                        if foundId not in addThermal:
                            if foundId not in presThermal:
                                addThermal.append(foundId)
                            if foundId not in addThermal:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addThermal.append(foundId)
                    reactNodes = _get_all_children(runNode, "react")
                    for reactNode in reactNodes:
                        sampleNode = _get_first_child(reactNode, "sample")
                        if sampleNode is not None:
                            foundId = sampleNode.attrib['id']
                            if foundId not in addSamples:
                                if foundId not in presSamples:
                                    addSamples.append(foundId)
                                if foundId not in addSamples:
                                    if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                        addSamples.append(foundId)
                        dataNodes = _get_all_children(reactNode, "data")
                        for dataNode in dataNodes:
                            lastNodes = _get_all_children(dataNode, "tar")
                            for lastNode in lastNodes:
                                foundId = lastNode.attrib['id']
                                if foundId not in addTargets:
                                    if foundId not in presTargets:
                                        addTargets.append(foundId)
                                    if foundId not in addTargets:
                                        if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                            addTargets.append(foundId)
                        partNodes = _get_all_children(reactNode, "partitions")
                        for partNode in partNodes:
                            dataNodes = _get_all_children(partNode, "data")
                            for dataNode in dataNodes:
                                lastNodes = _get_all_children(dataNode, "tar")
                                for lastNode in lastNodes:
                                    foundId = lastNode.attrib['id']
                                    if foundId not in addTargets:
                                        if foundId not in presTargets:
                                            addTargets.append(foundId)
                                        if foundId not in addTargets:
                                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                                addTargets.append(foundId)
        # Get the dependent of the dependent
        exp = _get_all_children(add_rd._node, "sample")
        for node in exp:
            if node.attrib['id'] in addSamples:
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)
                subNodes = _get_all_children(node, "type")
                for subNode in subNodes:
                    if 'targetId' in subNode.attrib:
                        foundId = subNode.attrib['targetId']
                        if foundId not in addTargets:
                            if foundId not in presTargets:
                                addTargets.append(foundId)
                            if foundId not in addTargets:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addTargets.append(foundId)
                subNodes = _get_all_children(node, "quantity")
                for subNode in subNodes:
                    if 'targetId' in subNode.attrib:
                        foundId = subNode.attrib['targetId']
                        if foundId not in addTargets:
                            if foundId not in presTargets:
                                addTargets.append(foundId)
                            if foundId not in addTargets:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addTargets.append(foundId)
                subNode = _get_first_child(node, "cdnaSynthesisMethod")
                if subNode is not None:
                    forId = _get_first_child(subNode, "thermalCyclingConditions")
                    if forId is not None:
                        foundId = forId.attrib['id']
                        if foundId not in addThermal:
                            if foundId not in presThermal:
                                addThermal.append(foundId)
                            if foundId not in addThermal:
                                if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                    addThermal.append(foundId)
        exp = _get_all_children(add_rd._node, "target")
        for node in exp:
            if node.attrib['id'] in addTargets:
                subNodes = _get_all_children(node, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addExps.append(foundId)
                subNode = _get_first_child(node, "dyeId")
                if subNode is not None:
                    foundId = subNode.attrib['id']
                    if foundId not in addDyes:
                        if foundId not in presDyes:
                            addDyes.append(foundId)
                        if foundId not in addDyes:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDyes.append(foundId)
        exp = _get_all_children(add_rd._node, "thermalCyclingConditions")
        for node in exp:
            if node.attrib['id'] in addThermal:
                subNodes = _get_all_children(node, "experimenter")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addExps:
                        if foundId not in presExps:
                            addExps.append(foundId)
                        if foundId not in addExps:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addExps.append(foundId)
                subNodes = _get_all_children(node, "documentation")
                for subNode in subNodes:
                    foundId = subNode.attrib['id']
                    if foundId not in addDocs:
                        if foundId not in presDocs:
                            addDocs.append(foundId)
                        if foundId not in addDocs:
                            if addMode in ["update-incl-dep", "all-incl-dep", "only-new-incl-dep"]:
                                addDocs.append(foundId)

        for docId in addExps:
            self.import_experimenter(add_rd.get_experimenter(byid=docId))
        for docId in addDocs:
            self.import_documentation(add_rd.get_documentation(byid=docId))
        for docId in addDyes:
            self.import_dye(add_rd.get_dye(byid=docId))
        for docId in addSamples:
            self.import_sample(add_rd, add_rd.get_sample(byid=docId), "no-dep")
        for docId in addTargets:
            self.import_target(add_rd, add_rd.get_target(byid=docId), "no-dep")
        for docId in addThermal:
            self.import_therm_cyc_cons(add_rd, add_rd.get_therm_cyc_cons(byid=docId), "no-dep")

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
        if key in ["description", "dyeChemistry", "dNTPs", "dyeConc"]:
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

        if key == "dyeChemistry":
            if value not in ["non-saturating DNA binding dye", "saturating DNA binding dye", "hybridization probe",
                             "hydrolysis probe", "labelled forward primer", "labelled reverse primer",
                             "DNA-zyme probe", ""]:
                raise RdmlError('Unknown or unsupported dye chemistry value "' + value + '".')

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        par = self._node.getparent()
        ver = par.get('version')
        if ver in ["1.3", "1.4"]:
            if key == "dyeChemistry":
                return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if ver in ["1.4"]:
            if key in ["dNTPs", "dyeConc"]:
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

        return ["id", "description", "dyeChemistry", "dNTPs", "dyeConc"]

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
        _add_first_child_to_dic(self._node, data, True, "dyeChemistry")
        _add_first_child_to_dic(self._node, data, True, "dNTPs")
        _add_first_child_to_dic(self._node, data, True, "dyeConc")
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
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key in ["interRunCalibrator", "doubleStranded", "calibratorSample"]:
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
        if ver in ["1.2", "1.3", "1.4"]:
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

        if key == "id":
            self.change_id(value, merge_with_id=False)
            return
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key in ["interRunCalibrator", "doubleStranded", "calibratorSample"]:
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
        if ver in ["1.2", "1.3", "1.4"]:
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
            return ["id", "description", "interRunCalibrator", "calibratorSample",
                    "cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                    "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions",
                    "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]
        return ["id", "description", "annotation", "interRunCalibrator", "doubleStranded", "calibratorSample",
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
                "doubleStranded", "quantity", "calibratorSample", "cdnaSynthesisMethod", "templateQuantity"]

    def types(self):
        """Returns a list of the types in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with type and id strings.
        """

        typesList = _get_all_children(self._node, "type")
        ret = []
        for node in typesList:
            data = {}
            data["type"] = node.text
            if "targetId" in node.attrib:
                data["targetId"] = node.attrib["targetId"]
            else:
                data["targetId"] = ""
            ret.append(data)
        return ret

    def new_type(self, type, targetId=None, newposition=None):
        """Creates a new type element.

        Args:
            self: The class self parameter.
            type: The "unkn", "ntc", "nac", "std", "ntp", "nrt", "pos" or "opt" type of sample
            targetId: The target linked to the type (makes sense in "pos" or "ntp" context)
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if type not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
            raise RdmlError('Unknown or unsupported sample type value "' + type + '".')
        new_node = et.Element("type")
        new_node.text = type
        par = self._node.getparent()
        ver = par.get('version')
        if ver in ["1.3", "1.4"]:
            if targetId is not None:
                if not targetId == "":
                    new_node.attrib["targetId"] = targetId
        place = _get_tag_pos(self._node, "type", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_type(self, type, oldposition, newposition=None, targetId=None):
        """Edits a type element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            type: The "unkn", "ntc", "nac", "std", "ntp", "nrt", "pos" or "opt" type of sample
            targetId: The target linked to the type (makes sense in "pos" or "ntp" context)

        Returns:
            Nothing, changes self.
        """

        if type not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
            raise RdmlError('Unknown or unsupported sample type value "' + type + '".')
        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a type.')

        pos = _get_tag_pos(self._node, "type", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "type", None, oldposition)
        ele.text = type
        par = self._node.getparent()
        ver = par.get('version')
        if "targetId" in ele.attrib:
            del ele.attrib["targetId"]
        if ver in ["1.3", "1.4"]:
            if targetId is not None:
                if not targetId == "":
                    ele.attrib["targetId"] = targetId
        self._node.insert(pos, ele)

    def move_type(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "type", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "type", None, oldposition)
        self._node.insert(pos, ele)

    def delete_type(self, byposition):
        """Deletes an type element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver in ["1.1", "1.2"]:
            ls = self.types()
            if len(ls) < 2:
                return

        elem = _get_first_child_by_pos_or_id(self._node, "type", None, byposition)
        self._node.remove(elem)

    def quantitys(self):
        """Returns a list of the types in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with type and id strings.
        """

        typesList = _get_all_children(self._node, "quantity")
        ret = []
        for node in typesList:
            data = {}
            _add_first_child_to_dic(node, data, False, "value")
            _add_first_child_to_dic(node, data, False, "unit")
            if "targetId" in node.attrib:
                data["targetId"] = node.attrib["targetId"]
            else:
                data["targetId"] = ""
            ret.append(data)
        return ret

    def new_quantity(self, value, unit, targetId=None, newposition=None):
        """Creates a new type element.

        Args:
            self: The class self parameter.
            value: The value of the quantity
            unit: The "cop", "fold", "dil", "ng", "nMol" or "other" type of quantity
            targetId: The target linked to the type (makes sense in "pos" or "ntp" context)
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if unit not in ["cop", "fold", "dil", "ng", "nMol", "other"]:
            raise RdmlError('Unknown or unsupported sample quantity unit "' + unit + '".')
        if value == "":
            raise RdmlError('Sample quantity value can not be empty.')

        par = self._node.getparent()
        ver = par.get('version')
        if ver in ["1.1", "1.2"]:
            ls = self.quantitys()
            if len(ls) != 0:
                return
        new_node = et.Element("quantity")
        _change_subelement(new_node, "value", ["value", "unit"], value, True, "float")
        _change_subelement(new_node, "unit", ["value", "unit"], unit, True, "string")
        if ver in ["1.3", "1.4"]:
            if targetId is not None:
                if not targetId == "":
                    new_node.attrib["targetId"] = targetId
        place = _get_tag_pos(self._node, "quantity", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_quantity(self, value, unit, oldposition, newposition=None, targetId=None):
        """Edits a type element.

        Args:
            self: The class self parameter.
            value: The value of the quantity
            unit: The "cop", "fold", "dil", "ng", "nMol" or "other" type of quantity
            oldposition: The old position of the element
            newposition: The new position of the element
            targetId: The target linked to the type (makes sense in "pos" or "ntp" context)

        Returns:
            Nothing, changes self.
        """

        if unit not in ["cop", "fold", "dil", "ng", "nMol", "other"]:
            raise RdmlError('Unknown or unsupported sample quantity unit "' + unit + '".')
        if value == "":
            raise RdmlError('Sample quantity value can not be empty.')
        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a quantity.')

        pos = _get_tag_pos(self._node, "quantity", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "quantity", None, oldposition)
        _change_subelement(ele, "value", ["value", "unit"], value, True, "float")
        _change_subelement(ele, "unit", ["value", "unit"], unit, True, "string")
        par = self._node.getparent()
        ver = par.get('version')
        if "targetId" in ele.attrib:
            del ele.attrib["targetId"]
        if ver in ["1.3", "1.4"]:
            if targetId is not None:
                if not targetId == "":
                    ele.attrib["targetId"] = targetId
        self._node.insert(pos, ele)

    def move_quantity(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "quantity", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "quantity", None, oldposition)
        self._node.insert(pos, ele)

    def delete_quantity(self, byposition):
        """Deletes an type element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "quantity", None, byposition)
        self._node.remove(elem)

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
        new_node = et.Element("xRef")
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
        new_node = et.Element("annotation")
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

    def edit_annotation_value(self, property, value=None):
        """Edits an annotation element.

        Args:
            self: The class self parameter.
            property: The property
            value: Its value

        Returns:
            Nothing, changes self.
        """

        par = self._node.getparent()
        ver = par.get('version')
        if ver == "1.1":
            return
        if property == "":
            raise RdmlError('A property is required to edit a annotation.')
        allAnnos = _get_all_children(self._node, "annotation")
        foundProp = False
        if value == None:
            foundProp = True
        if value == "":
            foundProp = True
        for node in allAnnos:
            curProp = _get_first_child_text(node, "property")
            if curProp == "":
               self._node.remove(node)
               continue
            if property == curProp:
                if value == None:
                    self._node.remove(node)
                    continue
                if value == "":
                    self._node.remove(node)
                    continue
                if foundProp == False:
                    _change_subelement(node, "value", ["property", "value"], value, True, "string")
                    foundProp = True
                else:
                    self._node.remove(node)
                    continue
        if foundProp == False:
           self.new_annotation(property=property, value=value, newposition=9999999999)

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
        if ver in ["1.2", "1.3", "1.4"]:
            data["annotations"] = self.annotations()
        data["types"] = self.types()
        _add_first_child_to_dic(self._node, data, True, "interRunCalibrator")
        if ver in ["1.4"]:
            _add_first_child_to_dic(self._node, data, True, "doubleStranded")
        data["quantitys"] = self.quantitys()
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
        if ver in ["1.2", "1.3", "1.4"]:
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
        if key in ["description", "amplificationEfficiencyMethod", "amplificationEfficiency",
                   "amplificationEfficiencySE", "meltingTemperature", "detectionLimit"]:
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
                   "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc",
                   "sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                   "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe1_oligoConc",
                   "sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                   "sequences_probe2_sequence", "sequences_probe2_oligoConc",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence", "sequences_amplicon_oligoConc"]:
            prim = _get_first_child(self._node, "sequences")
            if prim is None:
                return None
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc"]:
                sec = _get_first_child(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc"]:
                sec = _get_first_child(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                       "sequences_probe1_sequence", "sequences_probe1_oligoConc"]:
                sec = _get_first_child(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_probe2_sequence", "sequences_probe2_fivePrimeTag"]:
                sec = _get_first_child(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence", "sequences_amplicon_oligoConc"]:
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
            if key in ["sequences_forwardPrimer_oligoConc", "sequences_reversePrimer_oligoConc",
                       "sequences_probe1_oligoConc", "sequences_probe2_oligoConc",
                       "sequences_amplicon_oligoConc"]:
                return _get_first_child_text(sec, "oligoConc")
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
        if ver in ["1.2", "1.3", "1.4"]:
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

        par = self._node.getparent()
        ver = par.get('version')
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
        if ver in ["1.2", "1.3", "1.4"]:
            if key == "amplificationEfficiencySE":
                return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        if ver in ["1.3", "1.4"]:
            if key == "meltingTemperature":
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
                   "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc",
                   "sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                   "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe1_oligoConc",
                   "sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                   "sequences_probe2_sequence", "sequences_probe2_oligoConc",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence", "sequences_amplicon_oligoConc"]:
            prim = _get_or_create_subelement(self._node, "sequences", self.xmlkeys())
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc"]:
                sec = _get_or_create_subelement(prim, "forwardPrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc"]:
                sec = _get_or_create_subelement(prim, "reversePrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                       "sequences_probe1_sequence", "sequences_probe1_oligoConc"]:
                sec = _get_or_create_subelement(prim, "probe1",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_probe2_sequence", "sequences_probe2_oligoConc"]:
                sec = _get_or_create_subelement(prim, "probe2",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence", "sequences_amplicon_oligoConc"]:
                sec = _get_or_create_subelement(prim, "amplicon",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if sec is None:
                return None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_reversePrimer_threePrimeTag",
                       "sequences_probe1_threePrimeTag", "sequences_probe2_threePrimeTag",
                       "sequences_amplicon_threePrimeTag"]:
                _change_subelement(sec, "threePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence", "oligoConc"], value, True, "string")
            if key in ["sequences_forwardPrimer_fivePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_probe1_fivePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_amplicon_fivePrimeTag"]:
                _change_subelement(sec, "fivePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence", "oligoConc"], value, True, "string")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                _change_subelement(sec, "sequence",
                                   ["threePrimeTag", "fivePrimeTag", "sequence", "oligoConc"], value, True, "string")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                _change_subelement(sec, "sequence",
                                   ["threePrimeTag", "fivePrimeTag", "sequence", "oligoConc"], value, True, "string")
            if key in ["sequences_forwardPrimer_oligoConc", "sequences_reversePrimer_oligoConc",
                       "sequences_probe1_oligoConc", "sequences_probe2_oligoConc",
                       "sequences_amplicon_oligoConc"]:
                _change_subelement(sec, "oligoConc",
                                   ["threePrimeTag", "fivePrimeTag", "sequence", "oligoConc"], value, True, "string")

            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc"]:
                _remove_irrelevant_subelement(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc"]:
                _remove_irrelevant_subelement(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                       "sequences_probe1_sequence", "sequences_probe1_oligoConc"]:
                _remove_irrelevant_subelement(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_probe2_sequence", "sequences_probe2_oligoConc"]:
                _remove_irrelevant_subelement(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence", "sequences_amplicon_oligoConc"]:
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
        if ver in ["1.2", "1.3", "1.4"]:
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
            allExp = _get_all_children(par, "sample")
            for node in allExp:
                subNodes = _get_all_children(node, "type")
                for subNode in subNodes:
                    if "targetId" in subNode.attrib:
                        if subNode.attrib['targetId'] == oldValue:
                            subNode.attrib['targetId'] = value
                subNodes = _get_all_children(node, "quantity")
                for subNode in subNodes:
                    if "targetId" in subNode.attrib:
                        if subNode.attrib['targetId'] == oldValue:
                            subNode.attrib['targetId'] = value
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
                "amplificationEfficiencySE", "meltingTemperature", "detectionLimit", "dyeId",
                "sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                "sequences_forwardPrimer_sequence", "sequences_forwardPrimer_oligoConc",
                "sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                "sequences_reversePrimer_sequence", "sequences_reversePrimer_oligoConc",
                "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                "sequences_probe1_sequence", "sequences_probe1_oligoConc",
                "sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag",
                "sequences_probe2_sequence", "sequences_probe2_oligoConc",
                "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                "sequences_amplicon_sequence", "sequences_amplicon_oligoConc",
                "commercialAssay_company", "commercialAssay_orderNumber"]  # Also change in LinRegPCR save RDML

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "xRef", "type", "amplificationEfficiencyMethod",
                "amplificationEfficiency", "amplificationEfficiencySE", "meltingTemperature",
                "detectionLimit", "dyeId", "sequences", "commercialAssay"]

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
        new_node = et.Element("xRef")
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
        _add_first_child_to_dic(self._node, data, True, "amplificationEfficiencySE")
        _add_first_child_to_dic(self._node, data, True, "meltingTemperature")
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
                _add_first_child_to_dic(sec, sdic, True, "oligoConc")
                if len(sdic.keys()) != 0:
                    qdic["forwardPrimer"] = sdic
            sec = _get_first_child(elem, "reversePrimer")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                _add_first_child_to_dic(sec, sdic, True, "oligoConc")
                if len(sdic.keys()) != 0:
                    qdic["reversePrimer"] = sdic
            sec = _get_first_child(elem, "probe1")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                _add_first_child_to_dic(sec, sdic, True, "oligoConc")
                if len(sdic.keys()) != 0:
                    qdic["probe1"] = sdic
            sec = _get_first_child(elem, "probe2")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                _add_first_child_to_dic(sec, sdic, True, "oligoConc")
                if len(sdic.keys()) != 0:
                    qdic["probe2"] = sdic
            sec = _get_first_child(elem, "amplicon")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                _add_first_child_to_dic(sec, sdic, True, "oligoConc")
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
        new_node = et.Element("step")
        xml_temp_step = ["temperature", "duration", "temperatureChange", "durationChange", "measure", "ramp"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = et.SubElement(new_node, "temperature")
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
        new_node = et.Element("step")
        xml_temp_step = ["highTemperature", "lowTemperature", "duration", "temperatureChange",
                         "durationChange", "measure", "ramp"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = et.SubElement(new_node, "gradient")
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
        new_node = et.Element("step")
        xml_temp_step = ["goto", "repeat"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = et.SubElement(new_node, "loop")
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
        new_node = et.Element("step")
        xml_temp_step = ["temperature"]
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        subel = et.SubElement(new_node, "pause")
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
        new_node = et.Element("step")
        _add_new_subelement(new_node, "step", "nr", str(count + 1), False)
        et.SubElement(new_node, "lidOpen")
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
                if value not in ["", "real time", "meltcurve"]:
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

    def getreactjson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        allRuns = self.runs()
        runs = []
        runCount = 0
        for exp in allRuns:
            runs.append(exp.tojson())
            runs[runCount]["AllData"] = exp.getreactjson(curves=False)
            runCount += 1
        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["runs"] = runs

        ref_data = self.getExperimentData(refListOnly=True)
        data["used_references"] = ref_data["reference"]

        return data

    def getExperimentData(self, refListOnly=False):
        """Get the N0 data from an experiment.

        Args:
            self: The class self parameter.
            refListOnly: Return only a list of reference genes

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            references: A list of the references target ids
            N0: A dictionary with the results per target
        """

        res = {}
        res["reference"] = {}
        res["N0"] = {}
        res["Ncopy"] = {}
        res["Vol"] = {}
        tarType = {}
        refTar = {}

        allRuns = self.runs()
        pRoot = self._node.getparent()

        # Find all present Targets
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if "id" in target.attrib:
                tarId = target.attrib['id']
                tarType[tarId] = _get_first_child_text(target, "type")

        # Find all used Targets and N0
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                samId = _get_first_child(react, "sample")
                if samId is None:
                    continue
                if "id" not in samId.attrib:
                    continue
                sample = samId.attrib['id']
                volume = 20.0
                readVol = _get_first_child_text(react, "vol")
                if not readVol == "":
                    try:
                        readVol = float(readVol)
                    except ValueError:
                        pass
                    else:
                        if math.isfinite(readVol):
                            if readVol > 0.0:
                                volume =  readVol
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    target = tarId.attrib['id']
                    if sample not in res["Ncopy"]:
                        res["Ncopy"][sample] = {}
                    if sample not in res["Vol"]:
                        res["Vol"][sample] = {}
                    if target not in res["Ncopy"][sample]:
                        res["Ncopy"][sample][target] = []
                    if target not in res["Vol"][sample]:
                        res["Vol"][sample][target] = []
                    excluded = _get_first_child_text(react_data, "excl")
                    if excluded != "":
                        continue
                    if target not in tarType:
                        continue
                    if tarType[target] == "ref":
                        refTar[target] = 1
                    if refListOnly:
                        continue
                    nCopyVal = _get_first_child_text(react_data, "Ncopy")
                    if nCopyVal == "":
                        continue
                    try:
                        nCopyVal = float(nCopyVal)
                    except ValueError:
                        nCopyVal = -1.0
                    if not math.isfinite(nCopyVal):
                        nCopyVal = -1.0
                    corrFVal = _get_first_child_text(react_data, "corrF")
                    if corrFVal != "":
                        try:
                            corrFVal = float(corrFVal)
                        except ValueError:
                            corrFVal = 1.0
                        if math.isfinite(corrFVal):
                            nCopyVal *= corrFVal
                    corrPVal = _get_first_child_text(react_data, "corrP")
                    if corrPVal != "":
                        try:
                            corrPVal = float(corrPVal)
                        except ValueError:
                            corrPVal = 1.0
                        if math.isfinite(corrPVal):
                            if corrPVal != 0.0:
                                nCopyVal /= corrPVal
                    if math.isfinite(nCopyVal):
                        if nCopyVal > 0.0:
                            res["Ncopy"][sample][target].append(nCopyVal)
                            res["Vol"][sample][target].append(volume)

        res["reference"] = sorted(refTar, key=lambda key: refTar[key])
        return res

    def interRunCorr(self, overlapType="samples", selAnnotation="", updateRDML=False, calcCorrection=True):
        """Corrects inter run differences. Modifies the cq values and returns a json with additional data.

        Args:
            self: The class self parameter.
            overlapType: Base the overlap on "samples" or "annotation".
            selAnnotation: The annotation to use if overlapType == "annotation", else ignored.
            updateRDML: If true, update the RDML data with the calculated values.
            calcCorrection: If false, only the combined threshold and PCR efficiency is calculated

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            run: A list of the run ids
            target: A dictionary with the results per target
            plate: A dictionary with the results per plate
        """

        rootPar = self._node.getparent()
        dataVersion = rootPar.get('version')

        if dataVersion in ["1.0", "1.1", "1.2", "1.3"]:
            raise RdmlError('Error: InterRunCorrection requires at least RDML version 1.4. After upgrading version run LinRegPCR again.')

        res = {}
        if overlapType not in ["samples", "annotation"]:
            raise RdmlError('Error: Unknown overlap type.')
        if overlapType == "annotation":
            if selAnnotation == "":
                raise RdmlError('Error: Selection of annotation required.')
        res["runs"] = []
        res["target"] = {}
        res["plate"] = {}
        res["tsv"] = {}
        res["plate"]["corrP"] = []
        res["plate"]["matrix"] = []
        err = ""
        tarPara = {}
        samSel = {}
        allRuns = self.runs()

        # Get the sample infos
        pRoot = self._node.getparent()
        # Get Dic to look up sample type
        transSamTar = _sampleTypeToDics(pRoot)
        # Get Dic to look up annotation value by sample id
        if overlapType == "annotation":
            if selAnnotation != "":
                samples = _get_all_children(pRoot, "sample")
                for sample in samples:
                    if "id" in sample.attrib:
                        samId = sample.attrib['id']
                        xref = _get_all_children(sample, "annotation")
                        for node in xref:
                            anno = _get_first_child_text(node, "property")
                            if anno == selAnnotation:
                                val = _get_first_child_text(node, "value")
                                if val != "":
                                    samSel[samId] = val
        # Find all used Targets
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    tar = tarId.attrib['id']
                    res["target"][tar] = {}
                    res["target"][tar]["present"] = {}
                    res["target"][tar]["overlap"] = []

        sortTargets = sorted(list(res["target"].keys()))

        # Create Plate Matrix
        for row in range(0, len(allRuns)):
            res["plate"]["matrix"].append([])
            for col in range(0, len(allRuns)):
                res["plate"]["matrix"][row].append(-1.0)
                if row == col:
                    res["plate"]["matrix"][row][col] = -2.0
        res["plate"]["matrix"][0][0] = -2.0
        for tar in res["target"]:
            tarPara[tar] = {}
            for row in range(0, len(allRuns)):
                res["target"][tar]["overlap"].append([])
                for col in range(0, len(allRuns)):
                    res["target"][tar]["overlap"][row].append(0)
                    if row == col:
                        res["target"][tar]["overlap"][row][col] = -2.0

        # Collect Run - Target Information
        for pRunA in range(0, len(allRuns)):
            runA = allRuns[pRunA]
            res["runs"].append(runA['id'])
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    tar = tarId.attrib['id']
                    excluded = _get_first_child_text(react_data, "excl")
                    if excluded != "":
                        continue
                    res["target"][tar]["present"][pRunA] = True
 
        # Analyze the runs pair by pair
        if calcCorrection:
            for cRunA in range(0, len(allRuns)):
                condCount = {}
                for tar in res["target"]:
                    condCount[tar] = 0
                possible = {}
                runA = allRuns[cRunA]
                reacts = _get_all_children(runA._node, "react")
                for react in reacts:
                    samId = _get_first_child(react, "sample")
                    if samId is None:
                        continue
                    if "id" not in samId.attrib:
                        continue
                    sample = samId.attrib['id']
                    react_datas = _get_all_children(react, "data")
                    for react_data in react_datas:
                        tarId = _get_first_child(react_data, "tar")
                        if tarId is None:
                            continue
                        if "id" not in tarId.attrib:
                            continue
                        target = tarId.attrib['id']
                        if transSamTar[sample][target] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                            continue
                        excluded = _get_first_child_text(react_data, "excl")
                        if excluded != "":
                            continue
                        nCopyVal = _get_first_child_text(react_data, "Ncopy")
                        if nCopyVal == "":
                            continue
                        try:
                            nCopyVal = float(nCopyVal)
                        except ValueError:
                            continue
                        if not math.isfinite(nCopyVal):
                            continue
                        if nCopyVal <= 0.0:
                            continue
                        corrFVal = _get_first_child_text(react_data, "corrF")
                        if corrFVal != "":
                            try:
                                corrFVal = float(corrFVal)
                            except ValueError:
                                pass
                            if math.isfinite(corrFVal):
                                nCopyVal *= corrFVal

                        if target not in possible:
                            possible[target] = {}
                        if overlapType == "annotation":
                            if sample not in samSel:
                                continue
                            if samSel[sample] == "":
                                continue
                            translSamp = samSel[sample]
                            if translSamp not in possible[target]:
                                possible[target][translSamp] = []
                            possible[target][translSamp].append(nCopyVal)
                        else:
                            if sample not in possible[target]:
                                possible[target][sample] = []
                            possible[target][sample].append(nCopyVal)
                        condCount[target] += 1

                for cRunB in range(0, len(allRuns)):
                    bOverlap = {}
                    if cRunA == cRunB:
                        for tar in res["target"]:
                            res["target"][tar]["overlap"][cRunA][cRunA] = condCount[tar]
                        continue
                    if cRunA < cRunB:
                        continue
                    runB = allRuns[cRunB]
                    reacts = _get_all_children(runB._node, "react")
                    for react in reacts:
                        samId = _get_first_child(react, "sample")
                        if samId is None:
                            continue
                        if "id" not in samId.attrib:
                            continue
                        sample = samId.attrib['id']
                        react_datas = _get_all_children(react, "data")
                        for react_data in react_datas:
                            tarId = _get_first_child(react_data, "tar")
                            if tarId is None:
                                continue
                            if "id" not in tarId.attrib:
                                continue
                            target = tarId.attrib['id']
                            # Keep only overlapping values
                            if target not in possible:
                                continue
                            if overlapType == "annotation":
                                if sample not in samSel:
                                    continue
                                if samSel[sample] == "":
                                    continue
                                translSamp = samSel[sample]
                                if translSamp not in possible[target]:
                                    continue
                            else:
                                if sample not in possible[target]:
                                    continue
                            excluded = _get_first_child_text(react_data, "excl")
                            if excluded != "":
                                continue
                            nCopyVal = _get_first_child_text(react_data, "Ncopy")
                            if nCopyVal == "":
                                continue
                            try:
                                nCopyVal = float(nCopyVal)
                            except ValueError:
                                continue
                            if not math.isfinite(nCopyVal):
                                continue
                            if nCopyVal <= 0.0:
                                continue
                            corrFVal = _get_first_child_text(react_data, "corrF")
                            if corrFVal != "":
                                try:
                                    corrFVal = float(corrFVal)
                                except ValueError:
                                    pass
                                if math.isfinite(corrFVal):
                                    nCopyVal *= corrFVal

                            if target not in bOverlap:
                                bOverlap[target] = {}
                            if overlapType == "annotation":
                                if sample not in samSel:
                                    continue
                                if samSel[sample] == "":
                                    continue
                                translSamp = samSel[sample]
                                if translSamp not in bOverlap[target]:
                                    bOverlap[target][translSamp] = []
                                bOverlap[target][translSamp].append(nCopyVal)
                                res["target"][target]["overlap"][cRunA][cRunB] += 1
                                res["target"][target]["overlap"][cRunB][cRunA] += 1
                            else:
                                if sample not in bOverlap[target]:
                                    bOverlap[target][sample] = []
                                bOverlap[target][sample].append(nCopyVal)
                                res["target"][target]["overlap"][cRunA][cRunB] += 1
                                res["target"][target]["overlap"][cRunB][cRunA] += 1

                    plateSum = 0.0
                    plateNum = 0
                    for tar in bOverlap:
                        for sam in bOverlap[tar]:
                            for bVal in bOverlap[tar][sam]:
                                for aVal in possible[tar][sam]:
                                    plateSum += math.log(aVal/bVal)
                                    plateNum += 1
                        # npTarRes = np.ma.asarray(geoTarColl, dtype=np.float64)
                        # geoTarFac[tar] = scp.gmean(npTarRes)
                    if plateNum > 0:
                        plateMean = math.exp(plateSum / plateNum)
                    else:
                        plateMean = -1.0
                    if plateMean > 0.0:
                        res["plate"]["matrix"][cRunA][cRunB] = 1.0 / plateMean
                        res["plate"]["matrix"][cRunB][cRunA] = plateMean
                    else:
                        res["plate"]["matrix"][cRunA][cRunB] = -1.0
                        res["plate"]["matrix"][cRunB][cRunA] = -1.0

            # Set diagonal 1.0
            for sRun in range(0, len(allRuns)):
                maxVal = -3.0
                for curr in range(0, len(allRuns)):
                    maxVal = max(maxVal, res["plate"]["matrix"][sRun][curr])
                    maxVal = max(maxVal, res["plate"]["matrix"][curr][sRun])
                if maxVal > 0.0:
                    res["plate"]["matrix"][sRun][sRun] = 1.0

            for tar in sortTargets:
                for sRun in range(0, len(allRuns)):
                    if sRun in res["target"][tar]["present"] and res["target"][tar]["present"][sRun] is True:
                        pass
                    else:
                        for curr in range(0, len(allRuns)):
                            res["target"][tar]["overlap"][sRun][curr] = -10
                            res["target"][tar]["overlap"][curr][sRun] = -10

            # Handle a single run
            if len(allRuns) == 1:
                res["plate"]["matrix"][0][0] = 1.0

            # Fill matix gaps
            for appNr in range(0, 2):
                for mRow in range(0, len(allRuns)):
                    for mCol in range(0, len(allRuns)):
                        if -9.0 < res["plate"]["matrix"][mRow][mCol] < 0.0:
                            _pco_fixPlateMatix(res["plate"]["matrix"], mRow, mCol)
            for mRow in range(0, len(allRuns)):
                for mCol in range(0, len(allRuns)):
                    if -9.0 < res["plate"]["matrix"][mRow][mCol] < 0.0:
                        err = "Error Plate: Could not fix all matrix gaps."

            # Calc final correction factors
            for fRunA in range(0, len(allRuns)):
                linSum = 0.0
                linNum = 0
                for fRunB in range(0, len(allRuns)):
                    if res["plate"]["matrix"][fRunB][fRunA] > 0.0:
                        linSum += math.log(res["plate"]["matrix"][fRunB][fRunA])
                        linNum += 1
                if linNum > 0:
                    curRes = math.exp(linSum / linNum)
                    res["plate"]["corrP"].append(curRes)
                else:
                    res["plate"]["corrP"].append(-10.0)

        if err != "":
            res["error"] = err

        # Table creation
        runLine = ""
        for tRun in res["runs"]:
            runLine += "\t" + tRun

        res["tsv"]["run_correction_factors"] = runLine
        res["tsv"]["run_correction_factors"] += "\nCorrection Factor"
        for tCorr in res["plate"]["corrP"]:
            if tCorr > 0.0:
                res["tsv"]["run_correction_factors"] += "\t" + "{:.4f}".format(tCorr)
            else:
                res["tsv"]["run_correction_factors"] += "\tNot available"
        res["tsv"]["run_correction_factors"] += "\n"

        res["tsv"]["overlapping_conditions"] = ""
        for tar in sortTargets:
            res["tsv"]["overlapping_conditions"] += tar + runLine + "\n"
            for row in range(0, len(allRuns)):
                res["tsv"]["overlapping_conditions"] += res["runs"][row]
                for col in range(0, len(allRuns)):
                    overCount = res["target"][tar]["overlap"][row][col]
                    overOut = str(overCount)
                    if overCount < 0:
                        overOut = "Not available"
                    if overCount < -9:
                        overOut = "Not present"
                    res["tsv"]["overlapping_conditions"] += '\t' + overOut
                res["tsv"]["overlapping_conditions"] += '\n'
            for row in range(0, len(allRuns)):
                res["tsv"]["overlapping_conditions"] += '\t'
            res["tsv"]["overlapping_conditions"] += '\n'

        ##############################
        # write out the rdml results #
        ##############################
        if updateRDML is True:
            dataXMLelements = _getXMLDataType()
            for pRunA in range(0, len(allRuns)):
                runA = allRuns[pRunA]
                reacts = _get_all_children(runA._node, "react")
                for react in reacts:
                    react_datas = _get_all_children(react, "data")
                    for react_data in react_datas:
                        tarId = _get_first_child(react_data, "tar")
                        _change_subelement(react_data, "corrP", dataXMLelements, "-1.0", True, "string")
                        _change_subelement(react_data, "corrCq", dataXMLelements, "", True, "string")
                        corrPlate = -1.0
                        if tarId is None:
                            continue
                        if "id" not in tarId.attrib:
                            continue
                        tar = tarId.attrib['id']
                        corrPlate = res["plate"]["corrP"][pRunA]
                        goodVal = "{:.4f}".format(corrPlate)
                        _change_subelement(react_data, "corrP", dataXMLelements, goodVal, True, "string")
        return res

    def absoluteQuantification(self, method="reference", quantUnit="cop", estimate=True, reactionVolume=20.0, overlapType="samples", selAnnotation="", inclAnnotation=False):
        """Perform absolute quantification on the experiment.

        Args:
            self: The class self parameter.
            method: Base the overlap on "reference", "cq-guess" or "optical".
            quantUnit: The unit used for quantification "cop", "fold", "dil", "nMol", "ng" or "other"
            estimate: If true, missing targets are estimated.
            reactionVolume: The volume in microliter of a reation.
            overlapType: Base the overlap on "samples" or "annotation".
            selAnnotation: The annotation to use if overlapType == "annotation", else ignored.
            inclAnnotation: If true, all annotations are included in csv output

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
        """

        rootPar = self._node.getparent()
        dataVersion = rootPar.get('version')

        if dataVersion in ["1.0", "1.1", "1.2", "1.3"]:
            raise RdmlError('Error: AbsoluteQuantification requires at least RDML version 1.4. After upgrading version run LinRegPCR again.')

        if method not in ["reference", "cq-guess", "optical"]:
            raise RdmlError('Error: Unknown method used in absoluteQuantification.')
        if quantUnit not in ["cop", "fold", "dil", "nMol", "ng", "other"]:
            raise RdmlError('Error: Unknown quantUnit used in absoluteQuantification.')
        if overlapType not in ["samples", "annotation"]:
            raise RdmlError('Error: Unknown overlap type.')
        if overlapType == "annotation":
            if selAnnotation == "":
                raise RdmlError('Error: Selection of annotation required.')
        try:
            reactionVolume = float(reactionVolume)
        except ValueError:
            raise RdmlError('Error: Reaction volume must be a number.')
        if reactionVolume <= 0.0:
            raise RdmlError('Error: Reaction volume must be > 0.0.')

        copyThreshold_rot = _get_copies_threshold()
        mol = 6.02214076e23  # copies / Mole
        baseWeight = 660.0  # g / mol per base
        finalFluor = 1e12
        finalConc = -1.0
        opticalFactor = 1.0
        tresh = 0.0
        no_optical_calibrator = True
        nCopyData = self.getExperimentData()
        pRoot = self._node.getparent()
        transSamTar = _sampleTypeToDics(pRoot)

        res = {}
        tarType = {}
        samSelAnno = {}
        samAllAnnos = {}
        allAnnoKeys = {}
        overSelAnno = {}
        # corrNcopyFact corrects for standard and amplicon length
        #               but not for volume effects
        res["corrNcopyFact"] = {}
        res["absUnit"] = ""
        allRuns = self.runs()

        # Find the sample annotoation by sample id
        samples = _get_all_children(pRoot, "sample")
        for sample in samples:
            if "id" in sample.attrib:
                samId = sample.attrib['id']
                xref = _get_all_children(sample, "annotation")
                for node in xref:
                    anno = _get_first_child_text(node, "property")
                    if anno == "":
                        continue
                    val = _get_first_child_text(node, "value")
                    if val == "":
                        continue
                    if overlapType == "annotation":
                        if anno == selAnnotation:
                            samSelAnno[samId] = val
                    if anno not in allAnnoKeys:
                        allAnnoKeys[anno] = 1
                    if samId not in samAllAnnos:
                        samAllAnnos[samId] = {}
                    if anno not in samAllAnnos[samId]:
                        samAllAnnos[samId][anno] = val
        sortedAnnoKeys = sorted(allAnnoKeys.keys())

        # Check the annotation overlap
        if overlapType == "annotation":
            for samId in samAllAnnos:
                if selAnnotation not in samAllAnnos[samId]:
                    continue
                currSelVal = samAllAnnos[samId][selAnnotation]
                for selAnno in samAllAnnos[samId]:
                    if selAnno == selAnnotation:
                        continue
                    if selAnno not in overSelAnno:
                        overSelAnno[selAnno] = {}
                        overSelAnno[selAnno]["conf"] = False
                        overSelAnno[selAnno]["data"] = {}
                        overSelAnno[selAnno]["data"][currSelVal] = samAllAnnos[samId][selAnno]
                    else:
                        if currSelVal not in overSelAnno[selAnno]["data"]:
                            overSelAnno[selAnno]["data"][currSelVal] = samAllAnnos[samId][selAnno]
                        if overSelAnno[selAnno]["data"][currSelVal] != samAllAnnos[samId][selAnno]:
                            overSelAnno[selAnno]["conf"] = True

        # Find all target types tar: "toi"/"ref"
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if "id" in target.attrib:
                tarId = target.attrib['id']
                tarType[tarId] = _get_first_child_text(target, "type")

        # Find all used Targets and the geomean volume
        res["target"] = {}
        vol_sum = 0.0
        vol_num = 0
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                readVol = _get_first_child_text(react, "vol")
                if not readVol == "":
                    try:
                        readVol = float(readVol)
                    except ValueError:
                        pass
                    else:
                        if math.isfinite(readVol):
                            if readVol > 0.0:
                                vol_sum += readVol
                                vol_num += 1
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    tar = tarId.attrib['id']
                    res["target"][tar] = {}
        if vol_num > 0:
            meanVol = vol_sum / vol_num
            if meanVol > 0.0:
                reactionVolume = meanVol

        res["tsv"] = {}
        res["quantity"] = {}
        res["quantity"]["cop"] = {}
        res["quantity"]["cop"]["samples"] = {}
        res["quantity"]["cop"]["count"] = 0
        res["quantity"]["fold"] = {}
        res["quantity"]["fold"]["samples"] = {}
        res["quantity"]["fold"]["count"] = 0
        res["quantity"]["dil"] = {}
        res["quantity"]["dil"]["samples"] = {}
        res["quantity"]["dil"]["count"] = 0
        res["quantity"]["nMol"] = {}
        res["quantity"]["nMol"]["samples"] = {}
        res["quantity"]["nMol"]["count"] = 0
        res["quantity"]["ng"] = {}
        res["quantity"]["ng"]["samples"] = {}
        res["quantity"]["ng"]["count"] = 0
        res["quantity"]["other"] = {}
        res["quantity"]["other"]["samples"] = {}
        res["quantity"]["other"]["count"] = 0

        res["standard"] = {}
        for tar in res["target"]:
            res["corrNcopyFact"][tar] = -1.0
            # This are the used targets
            res["standard"][tar] = {}
            res["target"][tar]["ampliconLen"] = 100.0
            res["quantity"]["cop"]["samples"][tar] = {}
            res["quantity"]["fold"]["samples"][tar] = {}
            res["quantity"]["dil"]["samples"][tar] = {}
            res["quantity"]["nMol"]["samples"][tar] = {}
            res["quantity"]["ng"]["samples"][tar] = {}
            res["quantity"]["other"]["samples"][tar] = {}

        sortTargets = sorted(list(res["target"].keys()))

        # Collect the react data for samples used in the experiment sorted by targets
        stdSamples = {}
        optSamples = {}
        for pRun in range(0, len(allRuns)):
            runA = allRuns[pRun]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                samId = _get_first_child(react, "sample")
                if samId is None:
                    continue
                if "id" not in samId.attrib:
                    continue
                sampleID = samId.attrib['id']
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    target = tarId.attrib['id']
                    excluded = _get_first_child_text(react_data, "excl")
                    if excluded != "":
                        continue
                    if method == "optical":
                        if transSamTar[sampleID][target] != 'opt':
                            continue
                        else:
                            if sampleID not in optSamples:
                                optSamples[sampleID] = []
                            optSamples[sampleID].append(react_data)
                    else:
                        if transSamTar[sampleID][target] != 'std':
                            continue
                    if target not in stdSamples:
                        stdSamples[target] = {}
                    if sampleID not in stdSamples[target]:
                        stdSamples[target][sampleID] = []
                    stdSamples[target][sampleID].append(react_data)

        # Get the quantity information only for the used samples
        samples = _get_all_children(pRoot, "sample")
        for currTar in res["target"]:
            for sample in samples:
                if "id" in sample.attrib:
                    samId = sample.attrib['id']
                    if currTar not in stdSamples:
                        continue
                    if samId not in stdSamples[currTar]:
                        continue
                    subEles = _get_all_children(sample, "quantity")
                    for subNode in subEles:
                        if 'targetId' in subNode.attrib:
                            if subNode.attrib['targetId'] != currTar:
                                continue
                        qUnit = _get_first_child_text(subNode, "unit")
                        qValue = float(_get_first_child_text(subNode, "value"))
                        if qUnit == "dil":
                            if qValue != 0.0:
                                qValue = 1.0 / qValue
                        if qUnit in res["quantity"]:
                            if currTar in res["quantity"][qUnit]["samples"]:
                                if qValue not in res["quantity"][qUnit]["samples"][currTar]:
                                    res["quantity"][qUnit]["samples"][currTar][qValue] = {}
                                res["quantity"][qUnit]["samples"][currTar][qValue][samId] = 1
                                res["quantity"][qUnit]["count"] += 1

        optQuantity = {}
        for sample in samples:
            if "id" in sample.attrib:
                samId = sample.attrib['id']
                if samId not in optSamples:
                    continue
                subEles = _get_all_children(sample, "quantity")
                for subNode in subEles:
                    qUnit = _get_first_child_text(subNode, "unit")
                    if qUnit != "ng":
                        continue
                    qValue = float(_get_first_child_text(subNode, "value"))
                    if qValue not in optQuantity:
                        optQuantity[qValue] = {}
                    optQuantity[qValue][samId] = -1.0

        # Get the amplicon length only for the used targets
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if "id" in target.attrib:
                tarId = target.attrib['id']
                if tarId not in res["target"]:
                    continue
                seqEle = _get_first_child(target, "sequences")
                if seqEle is not None:
                    ampliconEle = _get_first_child(seqEle, "amplicon")
                    if ampliconEle is not None:
                        amplicon = _get_first_child_text(ampliconEle, "sequence")
                        if amplicon != "":
                            if len(amplicon) > 20:
                                res["target"][tarId]["ampliconLen"] = float(len(amplicon))

        if method == "reference":
            no_standard_found = True
            all_standards_found = True
            stdCurves = {}
            res["absUnit"] = quantUnit
            for qTar in res["target"]:
                currTarStd = False
                stdCurves[qTar] = {}
                for pRunA in range(0, len(allRuns)):
                    runA = allRuns[pRunA]
                    stdCurves[qTar][runA['id']] = {}
                    stdCurves[qTar][runA['id']]["x"] = []
                    stdCurves[qTar][runA['id']]["Ncopy"] = []
                if len(res["quantity"][quantUnit]["samples"][qTar]) > 0:
                    sortValues = sorted(list(res["quantity"][quantUnit]["samples"][qTar].keys()), reverse=True)
                    for qValue in res["quantity"][quantUnit]["samples"][qTar]:
                        if qValue <= 0.0:
                            continue
                        mean_sum = 0.0
                        mean_num = 0
                        vol_sum = 0.0
                        selSamples = list(res["quantity"][quantUnit]["samples"][qTar][qValue].keys())
                        for qSamp in selSamples:
                            for react_data in stdSamples[qTar][qSamp]:
                                corrFac = _get_first_child_text(react_data, "corrF")
                                calcCorr = 1.0
                                if not corrFac == "":
                                    try:
                                        calcCorr = float(corrFac)
                                    except ValueError:
                                        calcCorr = 1.0
                                    if not math.isfinite(calcCorr):
                                        calcCorr = 1.0
                                    if calcCorr > 1.0:
                                        calcCorr = 1.0
                                plateFac = _get_first_child_text(react_data, "corrP")
                                calcPlate = 1.0
                                if not plateFac == "":
                                    try:
                                        calcPlate = float(plateFac)
                                    except ValueError:
                                        calcPlate = 0.0
                                    if not math.isfinite(calcPlate):
                                        calcCorr = 0.0
                                    if calcPlate == 0.0:
                                        calcCorr = 0.0
                                    calcCorr = calcCorr / calcPlate
                                calcNcopy = _get_first_child_text(react_data, "Ncopy")
                                if calcCorr > 0.0001:
                                    if not calcNcopy == "":
                                        try:
                                            calcNcopy = float(calcNcopy)
                                        except ValueError:
                                            pass
                                        else:
                                            if math.isfinite(calcNcopy):
                                                finalNcopy = calcCorr * calcNcopy
                                                if finalNcopy > 0.0:
                                                    mean_sum += finalNcopy
                                                    mean_num += 1
                                                    volume = reactionVolume
                                                    react = react_data.getparent()
                                                    readVol = _get_first_child_text(react, "vol")
                                                    if not readVol == "":
                                                        try:
                                                            readVol = float(readVol)
                                                        except ValueError:
                                                            pass
                                                        else:
                                                            if math.isfinite(readVol):
                                                                if readVol > 0.0:
                                                                    volume = readVol
                                                    vol_sum += volume
                                                    run = react.getparent()
                                                    if 'id' in run.attrib:
                                                        stdCurves[qTar][run.attrib['id']]["x"].append(math.log10(qValue * volume))
                                                        stdCurves[qTar][run.attrib['id']]["Ncopy"].append(finalNcopy)
                        if qValue == sortValues[0]:
                            if mean_num > 0:
                                meanNcopy = mean_sum / mean_num
                                meanVol = vol_sum / mean_num
                                if meanVol > 0.0:
                                    reactionVolume = meanVol
                                res["corrNcopyFact"][qTar] = (sortValues[0] * 20.0) / meanNcopy  # (geoNcopy / sortValues[0]) * pow(1.9, 35.0)
                                no_standard_found = False
                                currTarStd = True
                    for pRunA in range(0, len(allRuns)):
                        for aPos  in range(0, len(stdCurves[qTar][run.attrib['id']]["Ncopy"])):
                            stdCurves[qTar][run.attrib['id']]["Ncopy"][aPos] = math.log2(copyThreshold_rot / (stdCurves[qTar][run.attrib['id']]["Ncopy"][aPos] * res["corrNcopyFact"][qTar]))
                if currTarStd == False:
                    all_standards_found = False
            if no_standard_found == False:
                if all_standards_found == False:
                    res["error"] = "Warning: Not all targets have standard samples. Consider to set estimate missing targets to yes."
            else:
                raise RdmlError("Error: No standard sample found. A sample type \"std\" must be selected and a quantity must be provided.")

            stdCurvesCsv = ""
            for tar in stdCurves:
                for runPos in range(0, len(allRuns)):
                    runEle = allRuns[runPos]
                    runID = runEle['id']
                    if runID in stdCurves[tar]:
                        if len(stdCurves[tar][runID]["x"]) > 1:
                            try:
                                curveLinReg = scp.stats.linregress(x=stdCurves[tar][runID]["x"], y=stdCurves[tar][runID]["Ncopy"])
                                dilPCREff = float(np.power(10.0, -1.0 / curveLinReg.slope))
                                rSquared = float(curveLinReg.rvalue) ** 2
                                dilFactor = 10 * float(np.power(2.0, -1.0 / np.log10(dilPCREff)))
                                stdCurvesCsv += tar + '\t'
                                stdCurvesCsv += runID + '\t'
                                stdCurvesCsv += "{:.4f}".format(rSquared) + '\t'
                                stdCurvesCsv += "{:.4f}".format(dilFactor) + '\n'
                            except ValueError:
                                stdCurvesCsv = ""

            if stdCurvesCsv != "":
                res["tsv"]["dilStandard"] = "Target\tRun\tDilution R^2\tDilution Factor\n"
                res["tsv"]["dilStandard"] += stdCurvesCsv

            if estimate is True:
                geoTar_sum = 0.0
                geoTar_num = 0
                for qTar in res["target"]:
                    if res["corrNcopyFact"][qTar] > 0.0:
                        geoTar_sum += math.log(res["corrNcopyFact"][qTar] * res["target"][tar]["ampliconLen"] / 100.0)
                        geoTar_num += 1
                if geoTar_num > 0:
                    geoTar = math.exp(geoTar_sum / geoTar_num)
                    for qTar in res["target"]:
                        if res["corrNcopyFact"][qTar] <= 0.0:
                            res["corrNcopyFact"][qTar] = geoTar * 100.0 / res["target"][tar]["ampliconLen"]

        if method == "cq-guess":
            res["absUnit"] = "cop"
            for tar in res["target"]:
                res["corrNcopyFact"][tar] = 100.0 / res["target"][tar]["ampliconLen"]

        if method == "optical":
            res["absUnit"] = "cop"
            concFluor = {}
            neg_sum = 0.0
            neg_num = 0
            tresh_sum = 0.0
            tresh_num = 0
            negFluor = 0.0
            no_optical_calibrator = True
            no_neg_calibrator = True
            if len(optQuantity) == 0:
                raise RdmlError('Error: No optical standard found.')
            for currConc in optQuantity:
                fluor_sum = 0.0
                fluor_num = 0
                for currSamp in optQuantity[currConc]:
                    for react_data in optSamples[currSamp]:
                        volume = reactionVolume
                        react = react_data.getparent()
                        readVol = _get_first_child_text(react, "vol")
                        if not readVol == "":
                            try:
                                readVol = float(readVol)
                            except ValueError:
                                pass
                            else:
                                if math.isfinite(readVol):
                                    if readVol > 0.0:
                                        volume = readVol
                        threshold = _get_first_child_text(react_data, "quantFluor")
                        if not threshold == "":
                            try:
                                threshold = float(threshold)
                            except ValueError:
                                continue
                            if math.isfinite(threshold):
                                if threshold > 0.0:
                                    tresh_sum += math.log(threshold)
                                    tresh_num += 1
                        adps = _get_all_children(react_data, "adp")
                        for adp in adps:
                            cyc = _get_first_child_text(adp, "cyc")
                            if not cyc == "":
                                try:
                                    cyc = int(float(cyc))
                                except ValueError:
                                    continue
                                if cyc < 0:
                                    continue
                            if 5 < cyc < 11:
                                curFlour = _get_first_child_text(adp, "fluor")
                                if not curFlour == "":
                                    try:
                                        curFlour = float(curFlour)
                                    except ValueError:
                                        continue
                                    if not math.isfinite(curFlour):
                                        continue
                                corrFlour = curFlour / volume
                                if currConc > 0.000001:
                                    fluor_sum += corrFlour
                                    fluor_num += 1
                                    no_optical_calibrator = False
                                else:
                                    neg_sum += corrFlour
                                    neg_num += 1
                                    no_neg_calibrator = False
                if fluor_num > 0:
                    concFluor[currConc] = fluor_sum / fluor_num
            if no_optical_calibrator == True:
                raise RdmlError("Error: Not optical calibrator sample found. A sample type \"opt\" must be selected and a quantity must be provided.")
            else:
                if no_neg_calibrator == True:
                    raise RdmlError("Error: Not negative optical calibrator sample found. A sample type \"opt\" must be selected and a 0.0 ng quantity must be provided.")

            if neg_num > 0:
                negFluor = neg_sum / neg_num
            if tresh_num > 0:
                tresh = math.exp(tresh_sum / tresh_num)
            if tresh <= 0.0:
                raise RdmlError('Error: No threshold found.')
            ngTreshold_rot = copyThreshold_rot * 1e9 * 100 * baseWeight / mol
            for currConc in concFluor:
                currFluor = (concFluor[currConc] / tresh) - (negFluor / tresh)
                if 0.05 + (negFluor / tresh) < currFluor < finalFluor:
                    finalFluor = currFluor
                    finalConc = currConc
            opticalFactor = finalConc / finalFluor / ngTreshold_rot
            for tar in res["target"]:
                res["corrNcopyFact"][tar] = opticalFactor * 100.0 / res["target"][tar]["ampliconLen"]
                # 1.0 / (float(res["threshold"]) / ((currConc * mol / (100 * baseWeight * 1e9)) / currFluor) * res["target"][tar]["ampliconLen"] / reactionVolume) / 100

        res["tsv"]["corrNcopyFact"] = 'Target\tQuant. Fact.\tAmplicon Length\t'
        res["tsv"]["corrNcopyFact"] += _niceQuantityTypeReact("cop") + ' at Threshold\t'
        res["tsv"]["corrNcopyFact"] += _niceQuantityTypeReact("ng") + ' at Threshold\t'
        res["tsv"]["corrNcopyFact"] += 'Reaction Volume\n'
        for tar in sortTargets:
            res["tsv"]["corrNcopyFact"] += tar + '\t'
            res["tsv"]["corrNcopyFact"] += "{:.6f}".format(res["corrNcopyFact"][tar]) + '\t'
            res["tsv"]["corrNcopyFact"] += str(res["target"][tar]["ampliconLen"]) + '\t'
            if res["corrNcopyFact"][tar] > 0.0:
                copiesCalc = copyThreshold_rot * res["corrNcopyFact"][tar] * reactionVolume / 20.0
                ngCalc = copiesCalc * 1e9 * res["target"][tar]["ampliconLen"] * baseWeight / mol
                res["tsv"]["corrNcopyFact"] += "{:.2e}".format(copiesCalc) + '\t'
                res["tsv"]["corrNcopyFact"] += "{:.2f}".format(ngCalc)
            else:
                res["tsv"]["corrNcopyFact"] += '\t'
            res["tsv"]["corrNcopyFact"] += '\t' + "{:.1f}".format(reactionVolume) + '\n'

        if method == "optical":
            csvStandard = ""
            if neg_num > 0:
                csvStandard += 'Threshold\t' + _niceQuantityTypeReact("ng") + '\t'
                csvStandard += "{:.2f}".format(tresh) + '\t\n'

            csvStandard += 'No DNA\t' + _niceQuantityTypeReact("ng") + '\t'
            csvStandard += "{:.2f}".format(negFluor * reactionVolume) + '\t\n'

            sortedConc = sorted(list(concFluor.keys()))
            ngTreshold_rot = copyThreshold_rot * 1e9 * 100 * baseWeight / mol
            for currConc in sortedConc:
                csvStandard += "{:.2f}".format(currConc * reactionVolume) + '\t' + _niceQuantityTypeReact("ng")
                resFluor = concFluor[currConc] / tresh - negFluor / tresh
                csvStandard += '\t' + "{:.2f}".format(resFluor * reactionVolume)
                corrFluor = resFluor * opticalFactor
                csvStandard += '\t' + "{:.2f}".format(corrFluor * reactionVolume)
                calcFluor =  currConc * reactionVolume / ngTreshold_rot
                csvStandard += '\t' + "{:.2f}".format(calcFluor) + '\n'

            if csvStandard != "":
                res["tsv"]["standard"] = "Concentration\tUnit\tRaw Fluor.\tCorrected Fluor.\tExpected Fluor.\n"
                res["tsv"]["standard"] += csvStandard
        else:
            csvStandard = ""
            for oUnit in ["cop", "fold", "dil", "nMol", "ng", "other"]:
                if res["absUnit"] != oUnit:
                    continue
                for oTar in res["quantity"][oUnit]["samples"]:
                    sortStdValues = sorted(list(res["quantity"][oUnit]["samples"][oTar].keys()), reverse=True)
                    for oValue in sortStdValues:
                        for oSample in res["quantity"][oUnit]["samples"][oTar][oValue]:
                            if oTar in stdSamples:
                                if oSample in stdSamples[oTar]:
                                    for react_data in stdSamples[oTar][oSample]:
                                        react = react_data.getparent()
                                        volume = reactionVolume
                                        readVol = _get_first_child_text(react, "vol")
                                        if not readVol == "":
                                            try:
                                                readVol = float(readVol)
                                            except ValueError:
                                                pass
                                            else:
                                                if math.isfinite(readVol):
                                                    if readVol > 0.0:
                                                        volume = readVol
                                        run = react.getparent()
                                        if 'id' in run.attrib:
                                            csvStandard += run.attrib['id']
                                        csvStandard += '\t'
                                        if 'id' in react.attrib:
                                            csvStandard += react.attrib['id']
                                        csvStandard += '\t' + oSample
                                        csvStandard += '\t' + oTar
                                        if oUnit == "cop":
                                            csvStandard += '\t' + "{:.2f}".format(oValue * volume) + '\t'
                                        else:
                                            csvStandard += '\t' + "{:.4e}".format(oValue * volume) + '\t'
                                        corrFac = _get_first_child_text(react_data, "corrF")
                                        calcCorr = 1.0
                                        if not corrFac == "":
                                            try:
                                                calcCorr = float(corrFac)
                                            except ValueError:
                                                calcCorr = 1.0
                                            if not math.isfinite(calcCorr):
                                                calcCorr = 1.0
                                            if calcCorr > 1.0:
                                                calcCorr = 1.0
                                        plateFac = _get_first_child_text(react_data, "corrP")
                                        calcPlate = 1.0
                                        if not plateFac == "":
                                            try:
                                                calcPlate = float(plateFac)
                                            except ValueError:
                                                calcPlate = 0.0
                                            if not math.isfinite(calcPlate):
                                                calcCorr = 0.0
                                            if calcPlate == 0.0:
                                                calcCorr = 0.0
                                            calcCorr = calcCorr / calcPlate
                                        calcNcopy = _get_first_child_text(react_data, "Ncopy")
                                        if calcCorr > 0.0001:
                                            if not calcNcopy == "":
                                                try:
                                                    calcNcopy = float(calcNcopy)
                                                except ValueError:
                                                    pass
                                                else:
                                                    if math.isfinite(calcNcopy):
                                                        finalNcopy = calcCorr * calcNcopy * volume / 20.0
                                                        calcQuant = finalNcopy * res["corrNcopyFact"][oTar]
                                                        if oUnit == "cop":
                                                            csvStandard += "{:.2f}".format(calcQuant)
                                                        else:
                                                            csvStandard += "{:.4e}".format(calcQuant)
                                        csvStandard += '\t' + _niceQuantityTypeReact(oUnit) + '\n'
            if csvStandard != "":
                res["tsv"]["standard"] = "Run\tReact\tSample\tTarget\tExpected\tCalculated\tUnit\n"
                res["tsv"]["standard"] += csvStandard

        # Mean the technical replicates
        res["tec_data"] = {}
        for sample in nCopyData["Ncopy"]:
            res["tec_data"][sample] = {}
            for target in nCopyData["Ncopy"][sample]:
                if transSamTar[sample][target] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                    continue
                res["tec_data"][sample][target] = {}
                res["tec_data"][sample][target]["sample_type"] = ""
                res["tec_data"][sample][target]["error"] = ""
                res["tec_data"][sample][target]["note"] = ""
                if sample in transSamTar:
                    if target in transSamTar[sample]:
                        res["tec_data"][sample][target]["sample_type"] = transSamTar[sample][target]
                res["tec_data"][sample][target]["target_type"] = ""
                if target in tarType:
                    res["tec_data"][sample][target]["target_type"] = tarType[target]
                res["tec_data"][sample][target]["n_tec_rep"] = len(nCopyData["Ncopy"][sample][target])
                res["tec_data"][sample][target]["raw_vals"] = []
                if len(nCopyData["Ncopy"][sample][target]) == 0:
                    res["tec_data"][sample][target]["error"] += "No Ncopy values;"
                    res["tec_data"][sample][target]["cop_mean"] = -1.0
                    res["tec_data"][sample][target]["cop_sd"] = -1.0
                    res["tec_data"][sample][target]["cop_cv"] = -1.0
                else:
                    for nPos in range(0, len(nCopyData["Ncopy"][sample][target])):
                        res["tec_data"][sample][target]["raw_vals"].append(nCopyData["Ncopy"][sample][target][nPos] * res["corrNcopyFact"][target] * nCopyData["Vol"][sample][target][nPos] / 20.0)
                    res["tec_data"][sample][target]["cop_mean"] = float(np.mean(res["tec_data"][sample][target]["raw_vals"]))
                    if len(nCopyData["Ncopy"][sample][target]) == 1:
                        res["tec_data"][sample][target]["cop_sd"] = 0.0
                        res["tec_data"][sample][target]["cop_cv"] = 0.0
                    else:
                        if res["tec_data"][sample][target]["cop_mean"] < 0.00:
                            res["tec_data"][sample][target]["error"] += "No Standard sample;"
                            res["tec_data"][sample][target]["cop_mean"] = -1.0
                            res["tec_data"][sample][target]["cop_sd"] = -1.0
                            res["tec_data"][sample][target]["cop_cv"] = -1.0
                        else:
                            res["tec_data"][sample][target]["cop_sd"] = float(np.std(res["tec_data"][sample][target]["raw_vals"], ddof=1))
                            calcCV = res["tec_data"][sample][target]["cop_sd"] / res["tec_data"][sample][target]["cop_mean"]
                            res["tec_data"][sample][target]["cop_cv"] = calcCV
                            if calcCV > 0.3:
                                res["tec_data"][sample][target]["note"] += "Tec. Rep. CV > 0.3;"

        if overlapType == "annotation":
            res["anno_data"] = {}
            res["anno_key"] = selAnnotation
            for sample in res["tec_data"]:
                for target in res["tec_data"][sample]:
                    if res["tec_data"][sample][target]["cop_mean"] > 0.0:
                        if sample in samSelAnno:
                            if target not in res["anno_data"]:
                                res["anno_data"][target] = {}
                            if samSelAnno[sample] not in res["anno_data"][target]:
                                res["anno_data"][target][samSelAnno[sample]] = {}
                            if "raw_vals" not in res["anno_data"][target][samSelAnno[sample]]:
                                res["anno_data"][target][samSelAnno[sample]]["raw_vals"] = []
                            res["anno_data"][target][samSelAnno[sample]]["raw_vals"].append(res["tec_data"][sample][target]["cop_mean"])
            for target in res["anno_data"]:
                for annoVal in res["anno_data"][target]:
                    annoCollVals = res["anno_data"][target][annoVal]["raw_vals"]
                    if len(annoCollVals) == 0:
                        res["anno_data"][target][annoVal]["mean"] = -1.0
                        res["anno_data"][target][annoVal]["sem"] = -1.0
                    else:
                        res["anno_data"][target][annoVal]["mean"] = float(np.mean(annoCollVals))
                        if len(annoCollVals) == 1:
                            res["anno_data"][target][annoVal]["sem"] = 0.0
                        else:
                            res["anno_data"][target][annoVal]["sem"] = float(scp.sem(annoCollVals))

        res["tsv"]["technical_data"] = "Sample\tSample Type\t"
        if inclAnnotation:
            for currAnno in sortedAnnoKeys:
                res["tsv"]["technical_data"] += currAnno + "\t"
        res["tsv"]["technical_data"] += "Target\tTarget Type\tError\tNote\t"
        res["tsv"]["technical_data"] += "n Tec. Rep.\tUnit\tMean\tSD\tCV\tIndividual Values\n"
        sortSam = sorted(res["tec_data"].keys())
        for sample in sortSam:
            sortTar = sorted(res["tec_data"][sample].keys())
            for target in sortTar:
                res["tsv"]["technical_data"] += sample + "\t"
                res["tsv"]["technical_data"] += res["tec_data"][sample][target]["sample_type"] + "\t"
                if inclAnnotation:
                    for currAnno in sortedAnnoKeys:
                        annoVal = ""
                        if sample in samAllAnnos:
                            if currAnno in samAllAnnos[sample]:
                                annoVal = samAllAnnos[sample][currAnno]
                        res["tsv"]["technical_data"] += annoVal + "\t"
                res["tsv"]["technical_data"] += target + "\t"
                res["tsv"]["technical_data"] += res["tec_data"][sample][target]["target_type"] + "\t"
                res["tsv"]["technical_data"] += res["tec_data"][sample][target]["error"] + "\t"
                res["tsv"]["technical_data"] += res["tec_data"][sample][target]["note"] + "\t"
                res["tsv"]["technical_data"] += str(res["tec_data"][sample][target]["n_tec_rep"]) + "\t"
                res["tsv"]["technical_data"] += res["absUnit"] + "\t"
                if res["absUnit"] == "cop":
                    res["tsv"]["technical_data"] += "{:.2f}".format(
                        res["tec_data"][sample][target]["cop_mean"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.2f}".format(
                        res["tec_data"][sample][target]["cop_sd"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.2f}".format(
                        res["tec_data"][sample][target]["cop_cv"]) + "\t"
                    if res["tec_data"][sample][target]["cop_mean"] > 0.0:
                        for indivVal in res["tec_data"][sample][target]["raw_vals"]:
                            res["tsv"]["technical_data"] += "{:.2f}".format(indivVal) + ";"
                else:
                    res["tsv"]["technical_data"] += "{:.6e}".format(
                        res["tec_data"][sample][target]["cop_mean"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.6e}".format(
                        res["tec_data"][sample][target]["cop_sd"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.6f}".format(
                        res["tec_data"][sample][target]["cop_cv"]) + "\t"
                    if res["tec_data"][sample][target]["cop_mean"] > 0.0:
                        for indivVal in res["tec_data"][sample][target]["raw_vals"]:
                            res["tsv"]["technical_data"] += "{:.4e}".format(indivVal) + ";"

                res["tsv"]["technical_data"] += "\n"

        if overlapType == "annotation":
            res["tsv"]["annotation_data"] = res["anno_key"] + "\t"
            if inclAnnotation:
                for currAnno in sortedAnnoKeys:
                    if currAnno == selAnnotation:
                        continue
                    if currAnno not in overSelAnno:
                        continue
                    if overSelAnno[currAnno]["conf"]:
                        continue
                    res["tsv"]["annotation_data"] += currAnno + "\t"
            res["tsv"]["annotation_data"] += "Target\tAbs. Expression\tSEM\tIndividual Values\n"
            sortTar = sorted(res["anno_data"].keys())
            for target in sortTar:
                sortAnno = sorted(res["anno_data"][target].keys())
                for annoVal in sortAnno:
                    res["tsv"]["annotation_data"] += annoVal + "\t"
                    if inclAnnotation:
                        for currAnno in sortedAnnoKeys:
                            if currAnno == selAnnotation:
                                continue
                            if currAnno not in overSelAnno:
                                continue
                            if overSelAnno[currAnno]["conf"]:
                                continue
                            annoData = ""
                            if currAnno in overSelAnno:
                                if annoVal in overSelAnno[currAnno]["data"]:
                                    annoData = overSelAnno[currAnno]["data"][annoVal]
                            res["tsv"]["annotation_data"] += annoData + "\t"
                    res["tsv"]["annotation_data"] += target + "\t"
                    res["tsv"]["annotation_data"] += "{:.4f}".format(res["anno_data"][target][annoVal]["mean"]) + "\t"
                    res["tsv"]["annotation_data"] += "{:.4f}".format(res["anno_data"][target][annoVal]["sem"]) + "\t"
                    for indivVal in res["anno_data"][target][annoVal]["raw_vals"]:
                        res["tsv"]["annotation_data"] += "{:.4f}".format(indivVal) + ";"
                    res["tsv"]["annotation_data"] += "\n"
        return res

    def genorm(self, selSamples="samples", selAnnotation="", selAnnoValue="", saveResultsCSV=False, saveResultsSVG=False, maxRef=-1):
        """Finds most stable reference genes. Returns a json with additional data.

        Args:
            self: The class self parameter.
            selSamples: Select "all" samples or samples with the "annotation" based on selAnnotation and selAnnoValue.
            selAnnotation: The annotation to use if overlapType == "annotation", else ignored.
            selAnnoValue: The value of the annotation to use if overlapType == "annotation", else ignored.
            saveResultsCSV: Save the results as tsv file
            saveResultsSVG: Save results as svg, this sets saveResultsCSV=True
            maxRef: Maximum number of reference genes allowed

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            run: A list of the run ids
            target: A dictionary with the results per target
            plate: A dictionary with the results per plate
        """

        rootPar = self._node.getparent()
        dataVersion = rootPar.get('version')

        if dataVersion in ["1.0", "1.1", "1.2", "1.3"]:
            raise RdmlError('Error: geNorm requires at least RDML version 1.4. After upgrading version run LinRegPCR again.')

        res = {}
        if selSamples not in ["all", "annotation"]:
            raise RdmlError('Error: Unknown sample selection.')
        if selSamples == "annotation":
            if selAnnotation == "":
                raise RdmlError('Error: Selection of annotation required.')
            if selAnnoValue == "":
                raise RdmlError('Error: Selection of annotation value required.')

        if saveResultsSVG:
            saveResultsCSV = True

        res["tsv"] = {}
        err = ""
        samSel = {}
        tarType = {}
        usedCond = {}
        usedTar = {}
        lookupCond = {}
        lookupTar = {}
        allRuns = self.runs()

        # Get the sample infos
        pRoot = self._node.getparent()
        transSamTar = _sampleTypeToDics(pRoot)
        if selSamples == "annotation":
            if selAnnotation != "":
                samples = _get_all_children(pRoot, "sample")
                for sample in samples:
                    if "id" in sample.attrib:
                        samId = sample.attrib['id']
                        xref = _get_all_children(sample, "annotation")
                        for node in xref:
                            anno = _get_first_child_text(node, "property")
                            if anno == selAnnotation:
                                val = _get_first_child_text(node, "value")
                                if val != "":
                                    samSel[samId] = val

        # Find all present Targets
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if "id" in target.attrib:
                tarId = target.attrib['id']
                tarType[tarId] = _get_first_child_text(target, "type")

        # Find all used Targets and Conditions
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                samId = _get_first_child(react, "sample")
                if samId is None:
                    continue
                if "id" not in samId.attrib:
                    continue
                sample = samId.attrib['id']
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    target = tarId.attrib['id']
                    if transSamTar[sample][target] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                        continue
                    excluded = _get_first_child_text(react_data, "excl")
                    if excluded != "":
                        continue
                    if tarType[target] == "ref":
                        usedTar[target] = 1
                        if selSamples == "annotation":
                            if sample not in samSel:
                                continue
                            if samSel[sample] == selAnnoValue:
                                usedCond[sample] = 1
                        else:
                            usedCond[sample] = 1

        res["reference"] = sorted(usedTar, key=lambda key: usedTar[key])
        res["conditions"] = sorted(usedCond, key=lambda key: usedCond[key])
        num = 0
        for cRef in res["reference"]:
            lookupTar[cRef] = num
            num += 1
        num = 0
        for cCon in res["conditions"]:
            lookupCond[cCon] = num
            num += 1

        if maxRef > 0:
            if len(res["reference"]) > maxRef:
                raise RdmlError('Server Error: geNorm can handle up to ' + str(maxRef) + ' reference genes.' + 
                                str(len(res["reference"])) + ' reference genes were given.')
        if len(res["reference"]) < 2:
            raise RdmlError('Error: geNorm requires at least two reference genes.')
        if len(res["conditions"]) < 2:
            raise RdmlError('Error: geNorm requires at least two conditions.')

        # Span the matrix
        gridSize = (len(res["conditions"]), len(res["reference"]))
        nCopy_sum = np.zeros(gridSize, dtype=np.float64)
        nCopy_num = np.zeros(gridSize, dtype=np.int64)
        nCopy_geo = np.zeros(gridSize, dtype=np.float64)

        # Fill the matrix
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                samId = _get_first_child(react, "sample")
                if samId is None:
                    continue
                if "id" not in samId.attrib:
                    continue
                sample = samId.attrib['id']
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    target = tarId.attrib['id']
                    if target not in res["reference"]:
                        continue
                    if sample not in res["conditions"]:
                        continue
                    nCopyVal = _get_first_child_text(react_data, "Ncopy")
                    if nCopyVal == "":
                        continue
                    try:
                        nCopyVal = float(nCopyVal)
                    except ValueError:
                        continue
                    if not math.isfinite(nCopyVal):
                        continue
                    if nCopyVal <= 0.0:
                        continue
                    corrFVal = _get_first_child_text(react_data, "corrF")
                    if corrFVal != "":
                        try:
                            corrFVal = float(corrFVal)
                        except ValueError:
                            pass
                        if math.isfinite(corrFVal):
                            nCopyVal *= corrFVal
                    corrPVal = _get_first_child_text(react_data, "corrP")
                    if corrPVal != "":
                        try:
                            corrPVal = float(corrPVal)
                        except ValueError:
                            pass
                        if math.isfinite(corrPVal):
                            if corrPVal != 0.0:
                                nCopyVal /= corrPVal
                    nCopy_sum[lookupCond[sample], lookupTar[target]] += np.log(nCopyVal)
                    nCopy_num[lookupCond[sample], lookupTar[target]] += 1
        with np.errstate(divide='ignore', invalid='ignore'):
            nCopy_geo = np.exp(nCopy_sum / nCopy_num)

        # Save the raw data
        res["reference"] = sorted(usedTar, key=lambda key: usedTar[key])
        if saveResultsCSV:
            res["tsv"]["nCopy_count"] = "\t"
            res["tsv"]["nCopy_values"] = "\t"
            for tar in res["reference"]:
                res["tsv"]["nCopy_count"] += tar + "\t"
                res["tsv"]["nCopy_values"] += tar + "\t"
            res["tsv"]["nCopy_count"] = re.sub(r"\t$", "\n", res["tsv"]["nCopy_count"])
            res["tsv"]["nCopy_values"] = re.sub(r"\t$", "\n", res["tsv"]["nCopy_values"])
            for row in range(0, len(res["conditions"])):
                res["tsv"]["nCopy_count"] += res["conditions"][row] + "\t"
                res["tsv"]["nCopy_values"] += res["conditions"][row] + "\t"
                for col in range(0, len(nCopy_num[row])):
                    res["tsv"]["nCopy_count"] += str(nCopy_num[row, col]) + "\t"
                    res["tsv"]["nCopy_values"] += "{:.4f}".format(nCopy_geo[row, col]) + "\t"
                res["tsv"]["nCopy_count"] = re.sub(r"\t$", "\n", res["tsv"]["nCopy_count"])
                res["tsv"]["nCopy_values"] = re.sub(r"\t$", "\n", res["tsv"]["nCopy_values"])

        # Remove empty columns and rows
        if np.amax(nCopy_num) == 0:
            res["error"] = "Error: No data to run geNorm+ on."
            return res
        columMax = np.amax(nCopy_num, axis=0)
        for col in range(len(columMax) -1, -1, -1):
            if columMax[col] == 0:
                err += "Removed target: " + res["reference"][col] + ";"
                del res["reference"][col]
                nCopy_num = np.delete(nCopy_num, col, axis=1)
                nCopy_geo = np.delete(nCopy_geo, col, axis=1)
        rowMax = np.amax(nCopy_num, axis=1)
        for row in range(len(rowMax) -1, -1, -1):
            if rowMax[row] == 0:
                err += "Removed condition: " + res["conditions"][row] + ";"
                del res["conditions"][row]
                nCopy_num = np.delete(nCopy_num, row, axis=0)
                nCopy_geo = np.delete(nCopy_geo, row, axis=0)

        # print(np.shape(nCopy_geo)[1])
        if np.shape(nCopy_geo)[1] < 2:
            raise RdmlError('Error: geNorm requires at least two reference genes.')
        if np.shape(nCopy_geo)[0] < 2:
            raise RdmlError('Error: geNorm requires at least two conditions.')

        # Calculate M factor
        mFactor = np.zeros(np.shape(nCopy_geo)[1], dtype=np.float64)
        for col in range(0, np.shape(nCopy_geo)[1]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                logRes = np.log2(nCopy_geo[:, col][:, None] / np.delete(nCopy_geo, col, axis=1))
                stdRes = np.nanstd(logRes, axis=0, ddof=1)
                mFactor[col] = np.nanmean(stdRes)
                if np.isnan(mFactor[col]):
                    raise RdmlError('Error: geNorm encountered not sufficient expression values for ' + res["reference"][col] + '.')
        #    print(res["reference"][col])
        #    print(stdRes)
        #    print(mFactor[col])

        # Calculate V factor
        vFactor = np.zeros(np.shape(nCopy_geo)[1] - 2, dtype=np.float64)
        mFactorInds = mFactor.argsort()
        sorted_mFactor = mFactor[mFactorInds[::1]]
        sorted_nCopy_geo = nCopy_geo[:, mFactorInds[::1]]
        sorted_Refs = np.array(res["reference"])[mFactorInds[::1]]

        res["result"] = ""
        res["m_targets"] = []
        res["m_values"] = []
        vValScore = -1
        vCount = 1
        for col in range(np.shape(nCopy_geo)[1] - 1, -1, -1):
            res["m_targets"].append(str(sorted_Refs[col]))
            res["m_values"].append(float(sorted_mFactor[col]))

        if saveResultsCSV:
            res["tsv"]["m_values"] = ""
            for tar in res["m_targets"]:
                res["tsv"]["m_values"] += tar + "\t"
            res["tsv"]["m_values"] = re.sub(r"\t$", "\n", res["tsv"]["m_values"])
            if np.shape(nCopy_geo)[1] > 2:
                for mVal in res["m_values"]:
                    res["tsv"]["m_values"] += "{:.6f}".format(mVal) + "\t"
                res["tsv"]["m_values"] = re.sub(r"\t$", "\n", res["tsv"]["m_values"])

        if np.shape(nCopy_geo)[1] > 2:
            res["v_labels"] = []
            res["v_values"] = []
            for col in range(1, np.shape(nCopy_geo)[1] - 1):
                res["v_labels"].append("v" + str(col+1) + "/" + str(col+2))
                sum_a = np.nansum(np.log(sorted_nCopy_geo[:, 0:col + 1]), axis=1)
                sum_b = np.nansum(np.log(sorted_nCopy_geo[:, 0:col + 2]), axis=1)
                num_a = np.count_nonzero(~np.isnan(np.log(sorted_nCopy_geo[:, 0:col + 1])), axis=1)
                num_b = np.count_nonzero(~np.isnan(np.log(sorted_nCopy_geo[:, 0:col + 2])), axis=1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    geoDiff = np.exp(sum_a / num_a) / np.exp(sum_b / num_b)
                    logDiff = np.log2(geoDiff)
                    vFactor[col - 1] = np.nanstd(logDiff, axis=0, ddof=1)
                    curr_vVal = float(np.nanstd(logDiff, axis=0, ddof=1))
                    res["v_values"].append(curr_vVal)
                    vCount += 1
                    if curr_vVal < 0.15:
                        if vValScore < 0:
                            vValScore = vCount

        if vValScore < 0:
            res["result"] = "No useful set of reference genes for normalization found (V >0.15). Start over by selecting new candidate reference genes."
        else:
            resString = "The combination of the genes "
            for genPos in range(0, vValScore):
                finPos = len(res["m_targets"]) - genPos - 1
                resString += res["m_targets"][finPos]
                if vValScore - 2 - genPos > 0:
                    resString += ", "
                if vValScore - 2 - genPos == 0:
                    resString += " and "
            resString += " is suitable for normalization."
            res["result"] = resString

        if saveResultsCSV:
            if np.shape(nCopy_geo)[1] > 2:
                res["tsv"]["v_values"] = ""
                for lab in res["v_labels"]:
                    res["tsv"]["v_values"] += lab + "\t"
                res["tsv"]["v_values"] = re.sub(r"\t$", "\n", res["tsv"]["v_values"])
                for vVal in res["v_values"]:
                    res["tsv"]["v_values"] += "{:.6f}".format(vVal) + "\t"
                res["tsv"]["v_values"] = re.sub(r"\t$", "\n", res["tsv"]["v_values"])

        if saveResultsSVG:
            res["svg"] = {}
            fig = plt_fig()
            axis = fig.add_subplot(1, 1, 1)
            axis.bar(res["m_targets"], res["m_values"])
            axis.tick_params(labelsize=16)
            axis.tick_params(axis="x", labelrotation=90)
            xLim = axis.get_xlim()
            yLim = axis.get_ylim()
            if yLim[1] < 0.8:
                axis.set_ylim(0.0, 0.8)
            else:
                axis.set_ylim(0.0, None)
            axis.plot(xLim, [0.5, 0.5], "k--")
            fig.tight_layout()
            with io.StringIO() as mFile:
                figSVG(fig).print_svg(mFile)
                mFile.seek(0)
                res["svg"]["m_values"] = mFile.read()

            if np.shape(nCopy_geo)[1] > 2:
                fig2 = plt_fig()
                axis2 = fig2.add_subplot(1, 1, 1)
                axis2.bar(res["v_labels"], res["v_values"])
                axis2.tick_params(labelsize=16)
                axis2.tick_params(axis="x", labelrotation=90)
                xLim2 = axis2.get_xlim()
                yLim2 = axis2.get_ylim()
                if yLim2[1] < 0.2:
                    axis2.set_ylim(0.0, 0.2)
                else:
                    axis2.set_ylim(0.0, None)
                axis2.plot(xLim2, [0.15, 0.15], "k--")
                fig2.tight_layout()
                with io.StringIO() as mFile2:
                    figSVG(fig2).print_svg(mFile2)
                    mFile2.seek(0)
                    res["svg"]["v_values"] = mFile2.read()

        if err != "":
            res["error"] = err

        return res

    def quantify(self, quantMethod="TD0", quantUnit="cop", overlapType="samples", selAnnotation="", statsParametric=False, statAlpha=0.05, inclAnnotation= False, selReferences=[], saveResultsCSV=False):
        """Calulates relative expression and returns a json with additional data.

        Args:
            self: The class self parameter.
            quantMethod: Quantify by "TD0" or "reference"
            quantUnit: The quantification unit "cop", "fold", "dil", "nMol", "ng" or "other"
            overlapType: Base the overlap on "samples" or "annotation".
            selAnnotation: The annotation to use if overlapType == "annotation", else ignored.
            statsParametric: True uses parametric tests, False non-parametric tests for statistics
            statAlpha: The alpha for group comparisons.
            inclAnnotation: If true, all annotations are included in csv output
            selReferences: The list of reference genes to correct for.
            saveResultsCSV: Save the results as tsv file.

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            run: A list of the run ids
            target: A dictionary with the results per target
            plate: A dictionary with the results per plate
        """

        if quantMethod not in ["TD0", "reference"]:
            quantMethod = "TD0"

        if quantUnit not in ["cop", "fold", "dil", "nMol", "ng", "other"]:
            quantUnit = "cop"

        rootPar = self._node.getparent()
        dataVersion = rootPar.get('version')

        if dataVersion in ["1.0", "1.1", "1.2", "1.3"]:
            raise RdmlError('Error: RelativeQuantification requires at least RDML version 1.4. After upgrading version run LinRegPCR again.')

       # if len(selReferences) < 1:
      #      raise RdmlError('Error: RelativeQuantification requires at least one reference gene.')

        res = {}
        tarType = {}
        samSelAnno = {}
        samAllAnnos = {}
        allAnnoKeys = {}
        overSelAnno = {}
        allRuns = self.runs()
        if overlapType not in ["samples", "annotation"]:
            raise RdmlError('Error: Unknown overlap type.')
        if overlapType == "annotation":
            if selAnnotation == "":
                raise RdmlError('Error: Selection of annotation required.')

        res["tec_data"] = {}
        res["ref_data"] = {}
        res["rel_data"] = {}
        res["tsv"] = {}
        err = ""

        nCopyData = self.getExperimentData()
        pRoot = self._node.getparent()
        transSamTar = _sampleTypeToDics(pRoot)

        reactionVolume = 10.0

        # Find the sample annotoation
        samples = _get_all_children(pRoot, "sample")
        for sample in samples:
            if "id" in sample.attrib:
                samId = sample.attrib['id']
                xref = _get_all_children(sample, "annotation")
                for node in xref:
                    anno = _get_first_child_text(node, "property")
                    if anno == "":
                        continue
                    val = _get_first_child_text(node, "value")
                    if val == "":
                        continue
                    if overlapType == "annotation":
                        if anno == selAnnotation:
                            samSelAnno[samId] = val
                    if anno not in allAnnoKeys:
                        allAnnoKeys[anno] = 1
                    if samId not in samAllAnnos:
                        samAllAnnos[samId] = {}
                    if anno not in samAllAnnos[samId]:
                        samAllAnnos[samId][anno] = val
        sortedAnnoKeys = sorted(allAnnoKeys.keys())

        # Check the annotation overlap
        if overlapType == "annotation":
            for samId in samAllAnnos:
                if selAnnotation not in samAllAnnos[samId]:
                    continue
                currSelVal = samAllAnnos[samId][selAnnotation]
                for selAnno in samAllAnnos[samId]:
                    if selAnno == selAnnotation:
                        continue
                    if selAnno not in overSelAnno:
                        overSelAnno[selAnno] = {}
                        overSelAnno[selAnno]["conf"] = False
                        overSelAnno[selAnno]["data"] = {}
                        overSelAnno[selAnno]["data"][currSelVal] = samAllAnnos[samId][selAnno]
                    else:
                        if currSelVal not in overSelAnno[selAnno]["data"]:
                            overSelAnno[selAnno]["data"][currSelVal] = samAllAnnos[samId][selAnno]
                        if overSelAnno[selAnno]["data"][currSelVal] != samAllAnnos[samId][selAnno]:
                            overSelAnno[selAnno]["conf"] = True

        # Find all target types
        scalingFact = {}
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if "id" in target.attrib:
                tarId = target.attrib['id']
                tarType[tarId] = _get_first_child_text(target, "type")
                scalingFact[tarId] = 1.0

        stdCurves = {}
        vol_sum = 0.0
        vol_num = 0
        for tRunA in range(0, len(allRuns)):
            runA = allRuns[tRunA]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                readVol = _get_first_child_text(react, "vol")
                if not readVol == "":
                    try:
                        readVol = float(readVol)
                    except ValueError:
                        pass
                    else:
                        if math.isfinite(readVol):
                            if readVol > 0.0:
                                vol_sum += readVol
                                vol_num += 1
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    tar = tarId.attrib['id']
                    stdCurves[tar] = {}
        if vol_num > 0:
            meanVol = vol_sum / vol_num
            if meanVol > 0.0:
                reactionVolume = meanVol

        # Analyze dilution curves
        no_standard_found = True
        all_standards_found = True
        sampleToQuantity = {}
        res["absUnit"] = quantUnit
        std_data = {}

        # Collect conc of standard samples with the used unit
        samples = _get_all_children(pRoot, "sample")
        for currTar in stdCurves:
            for sample in samples:
                if "id" in sample.attrib:
                    samId = sample.attrib['id']
                    subEles = _get_all_children(sample, "quantity")
                    for subNode in subEles:
                        if 'targetId' in subNode.attrib:
                            if subNode.attrib['targetId'] != currTar:
                                continue
                        qUnit = _get_first_child_text(subNode, "unit")
                        if qUnit != quantUnit:
                            continue
                        qValue = _get_first_child_text(subNode, "value")
                        try:
                            qValueFloat = float(qValue)
                        except ValueError:
                            pass
                        else:
                            if math.isfinite(qValueFloat):
                                if qUnit == "dil":
                                    if qValueFloat != 0.0:
                                        qValueFloat = 1.0 / qValueFloat
                                if currTar not in sampleToQuantity:
                                    sampleToQuantity[currTar] = {}
                                sampleToQuantity[currTar][samId] = qValueFloat
                                std_data[currTar] = {}
                                std_data[currTar]["raw"] = {}
                                std_data[currTar]["Cq"] = {}

        # Collect all Data
        pcrEff_sum = {}
        pcrEff_num = {}
        for pRun in range(0, len(allRuns)):
            runA = allRuns[pRun]
            reacts = _get_all_children(runA._node, "react")
            for react in reacts:
                samId = _get_first_child(react, "sample")
                if samId is None:
                    continue
                if "id" not in samId.attrib:
                    continue
                sampleID = samId.attrib['id']
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    tarId = _get_first_child(react_data, "tar")
                    if tarId is None:
                        continue
                    if "id" not in tarId.attrib:
                        continue
                    target = tarId.attrib['id']
                    ampEff = _get_first_child_text(react_data, "ampEff")
                    if not ampEff == "":
                        try:
                            ampEff = float(ampEff)
                        except ValueError:
                            pass
                        if math.isfinite(ampEff):
                            if target not in pcrEff_sum:
                                pcrEff_sum[target] = 0.0
                                pcrEff_num[target] = 0
                            pcrEff_sum[target] += ampEff
                            pcrEff_num[target] += 1
                    if target not in sampleToQuantity:
                        continue
                    if sampleID not in sampleToQuantity[target]:
                        continue
                    excluded = _get_first_child_text(react_data, "excl")
                    if excluded != "":
                        continue
                    corrFac = _get_first_child_text(react_data, "corrF")
                    calcCorr = 1.0
                    if not corrFac == "":
                        try:
                            calcCorr = float(corrFac)
                        except ValueError:
                            calcCorr = 1.0
                        if not math.isfinite(calcCorr):
                            calcCorr = 1.0
                        if calcCorr > 1.0:
                            calcCorr = 1.0
                    plateFac = _get_first_child_text(react_data, "corrP")
                    calcPlate = 1.0
                    if not plateFac == "":
                        try:
                            calcPlate = float(plateFac)
                        except ValueError:
                            calcPlate = 0.0
                        if not math.isfinite(calcPlate):
                            calcCorr = 0.0
                        if calcPlate == 0.0:
                            calcCorr = 0.0
                        calcCorr = calcCorr / calcPlate
                    calcNcopy = _get_first_child_text(react_data, "Ncopy")
                    cqStr = _get_first_child_text(react_data, "cq")
                    if calcCorr > 0.0001:
                        if not calcNcopy == "":
                            try:
                                calcNcopy = float(calcNcopy)
                            except ValueError:
                                pass
                            else:
                                if math.isfinite(calcNcopy):
                                    if not cqStr == "":
                                        try:
                                            cqFloat = float(cqStr)
                                        except ValueError:
                                            pass
                                        else:
                                            if math.isfinite(cqFloat):
                                                finalNcopy = calcCorr * calcNcopy
                                                if finalNcopy > 0.0:
                                                    volume = reactionVolume
                                                    react = react_data.getparent()
                                                    readVol = _get_first_child_text(react, "vol")
                                                    if not readVol == "":
                                                        try:
                                                            readVol = float(readVol)
                                                        except ValueError:
                                                            pass
                                                        else:
                                                            if math.isfinite(readVol):
                                                                if readVol > 0.0:
                                                                    volume = readVol
                                                    if quantUnit in ["cop", "nMol", "ng"]:
                                                        standardFix = str(sampleToQuantity[target][sampleID] * volume)
                                                    else:
                                                        standardFix = str(sampleToQuantity[target][sampleID])
                                                    # The data are as data{target}{"raw"}{"standard conc as string"}["float, "float"]
                                                    if standardFix not in std_data[target]["raw"]:
                                                        std_data[target]["raw"][standardFix] = []
                                                    if standardFix not in std_data[target]["Cq"]:
                                                        std_data[target]["Cq"][standardFix] = []
                                                    std_data[target]["raw"][standardFix].append(calcNcopy)
                                                    std_data[target]["Cq"][standardFix].append(cqFloat)
        std_data = standardCurveStats(std_data, False)

        resTars = sorted(list(std_data.keys()))
        res["dil_std"] = {}
        for tar in resTars:
            # Update scaling factor
            if quantMethod == "reference":
                if std_data[tar]["mean_max"] > 0.0 :
                    scalingFact[tar] = std_data[tar]["expected_max"] / std_data[tar]["mean_max"]

            if std_data[tar]["expected_min"] < 0.0:
                continue
            res["dil_std"][tar] = {}
            res["dil_std"][tar]["curve_pcr_eff"] = -1.0
            if tar in pcrEff_sum:
                res["dil_std"][tar]["curve_pcr_eff"] = pcrEff_sum[tar] / pcrEff_num[tar]
            res["dil_std"][tar]["dilution_pcr_eff"] = std_data[tar]["cq_eff"]
            res["dil_std"][tar]["mean_min"] = std_data[tar]["mean_min"]
            res["dil_std"][tar]["mean_max"] = std_data[tar]["mean_max"]
            res["dil_std"][tar]["mean_ratio"] = std_data[tar]["bias_as_ratio"]
            res["dil_std"][tar]["expected_min"] = std_data[tar]["expected_min"]
            res["dil_std"][tar]["expected_max"] = std_data[tar]["expected_max"]
            res["dil_std"][tar]["expected_ratio"] = std_data[tar]["expected_bias_as_ratio"]
            res["dil_std"][tar]["slope_bias"] = std_data[tar]["slope_bias"]
            res["dil_std"][tar]["correlation_R"] = std_data[tar]["correlation_R"]
            res["dil_std"][tar]["linearity"] = std_data[tar]["MSderivation"]
            res["dil_std"][tar]["reproducibility"] = std_data[tar]["MSwithin"]
            res["dil_std"][tar]["ms_within_linregpcr"] = std_data[tar]["cq_ms_within_linregpcr"]
            res["dil_std"][tar]["ratio"] = std_data[tar]["cq_ms_ratio"]
            res["dil_std"][tar]["detectable_diff"] = std_data[tar]["dif_detectable_diff"]

        # Scale the samples if requested
        if quantMethod == "reference":
            for sample in nCopyData["Ncopy"]:
                res["tec_data"][sample] = {}
                for target in nCopyData["Ncopy"][sample]:
                    if target in scalingFact:
                        for pos in range(0, len(nCopyData["Ncopy"][sample][target])):
                            nCopyData["Ncopy"][sample][target][pos] *= float(scalingFact[target])

        # Mean the technical replicates
        stats_tec_replicates_n = 0
        stats_tec_replicates_count = 0
        stats_tec_replicates_sum = 0.0
        stats_tec_replicates_geo = 0.0

        for sample in nCopyData["Ncopy"]:
            res["tec_data"][sample] = {}
            for target in nCopyData["Ncopy"][sample]:
                if transSamTar[sample][target] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                    continue
                res["tec_data"][sample][target] = {}
                res["tec_data"][sample][target]["sample_type"] = ""
                res["tec_data"][sample][target]["error"] = ""
                res["tec_data"][sample][target]["note"] = ""
                if sample in transSamTar:
                    if target in transSamTar[sample]:
                        res["tec_data"][sample][target]["sample_type"] = transSamTar[sample][target]
                res["tec_data"][sample][target]["target_type"] = ""
                if target in tarType:
                    res["tec_data"][sample][target]["target_type"] = tarType[target]
                res["tec_data"][sample][target]["n_tec_rep"] = len(nCopyData["Ncopy"][sample][target])
                res["tec_data"][sample][target]["raw_vals"] = nCopyData["Ncopy"][sample][target]
                if len(nCopyData["Ncopy"][sample][target]) == 0:
                    res["tec_data"][sample][target]["error"] += "No Ncopy values;"
                    res["tec_data"][sample][target]["Ncopy_mean"] = -1.0
                    res["tec_data"][sample][target]["Ncopy_sd"] = -1.0
                    res["tec_data"][sample][target]["Ncopy_cv"] = -1.0
                else:
                    res["tec_data"][sample][target]["Ncopy_mean"] = float(np.mean(nCopyData["Ncopy"][sample][target]))
                    stats_tec_replicates_min = float(np.min(nCopyData["Ncopy"][sample][target]))
                    if len(nCopyData["Ncopy"][sample][target]) == 1:
                        res["tec_data"][sample][target]["Ncopy_sd"] = 0.0
                        res["tec_data"][sample][target]["Ncopy_cv"] = 0.0
                    else:
                        res["tec_data"][sample][target]["Ncopy_sd"] = float(np.std(nCopyData["Ncopy"][sample][target], ddof=1))
                        calcCV = res["tec_data"][sample][target]["Ncopy_sd"] / res["tec_data"][sample][target]["Ncopy_mean"]
                        res["tec_data"][sample][target]["Ncopy_cv"] = calcCV
                        if stats_tec_replicates_min > 50:
                            if calcCV > 0.0:
                                stats_tec_replicates_n += 1
                                stats_tec_replicates_count += len(nCopyData["Ncopy"][sample][target])
                                stats_tec_replicates_sum += calcCV
                                stats_tec_replicates_geo += math.log(calcCV)
                        if calcCV > 0.3:
                            res["tec_data"][sample][target]["note"] += "Tec. Rep. CV > 0.3;"

        res["stat_tec"] = {}
        res["stat_tec"]["tec_replicates_n"] = stats_tec_replicates_n
        res["stat_tec"]["tec_replicates_aver"] = float(stats_tec_replicates_count) / float(stats_tec_replicates_n)
        res["stat_tec"]["tec_replicates_mean_SV"] = 0.0
        res["stat_tec"]["tec_replicates_geomean_SV"] = 0.0
        if stats_tec_replicates_n > 0:
            res["stat_tec"]["tec_replicates_mean_SV"] = stats_tec_replicates_sum / float(stats_tec_replicates_n)
            res["stat_tec"]["tec_replicates_geomean_SV"] = math.exp(stats_tec_replicates_geo / float(stats_tec_replicates_n))

        # TODO: remove
        print(res["stat_tec"]["tec_replicates_n"])
        print(res["stat_tec"]["tec_replicates_aver"])
        print(res["stat_tec"]["tec_replicates_mean_SV"])
        print(res["stat_tec"]["tec_replicates_geomean_SV"])

        if len(selReferences) > 0:
            # Geomean the reference genes
            refGeoSum = 0.0
            refGeoNum = 0
            for sample in nCopyData["Ncopy"]:
                keepSample = False
                for target in selReferences:
                    if transSamTar[sample][target] not in ["ntc", "nac", "ntp", "nrt", "opt"]:
                        keepSample = True
                if not keepSample:
                    continue
                res["ref_data"][sample] = {}
                res["ref_data"][sample]["Ncopy_sum"] = 0.0
                res["ref_data"][sample]["Ncopy_num"] = 0
                res["ref_data"][sample]["Ncopy_gem"] = -1.0
                res["ref_data"][sample]["ref_missing"] = False
                res["ref_data"][sample]["raw_vals"] = []
                for target in selReferences:
                    if sample not in res["tec_data"]:
                        continue
                    if target not in res["tec_data"][sample]:
                        res["ref_data"][sample]["ref_missing"] = True
                        continue
                    if res["tec_data"][sample][target]["n_tec_rep"] == 0:
                        res["ref_data"][sample]["ref_missing"] = True
                        continue
                    res["ref_data"][sample]["raw_vals"].append(res["tec_data"][sample][target]["Ncopy_mean"])
                    res["ref_data"][sample]["Ncopy_sum"] += math.log(res["tec_data"][sample][target]["Ncopy_mean"])
                    refGeoSum += math.log(res["tec_data"][sample][target]["Ncopy_mean"])
                    res["ref_data"][sample]["Ncopy_num"] += 1
                    refGeoNum += 1
                if res["ref_data"][sample]["Ncopy_num"] > 0:
                    geoMean = math.exp(res["ref_data"][sample]["Ncopy_sum"] / res["ref_data"][sample]["Ncopy_num"])
                    res["ref_data"][sample]["Ncopy_gem"] = geoMean
            res["ref_geoAll"] = 1.0
            if refGeoNum > 0:
                refGeoAll = math.exp(refGeoSum / refGeoNum)
                res["ref_geoAll"] = refGeoAll

            # Calculate relative gene expression
            for sample in nCopyData["Ncopy"]:
                res["rel_data"][sample] = {}
                for target in nCopyData["Ncopy"][sample]:
                    if sample not in transSamTar:
                        continue
                    if target not in transSamTar[sample]:
                        continue
                    if transSamTar[sample][target] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                        continue
                    if sample not in res["tec_data"]:
                        continue
                    if target not in res["tec_data"][sample]:
                        continue
                    if "target_type" not in res["tec_data"][sample][target]:
                        continue
                    if res["tec_data"][sample][target]["target_type"] == "ref":
                        continue
                    res["rel_data"][sample][target] = {}
                    res["rel_data"][sample][target]["rel_expression"] = -1.0
                    res["rel_data"][sample][target]["ref_missing"] = res["ref_data"][sample]["ref_missing"]
                    if not res["rel_data"][sample][target]["ref_missing"]:
                        if sample not in res["ref_data"]:
                            continue
                        relEx = -1.0
                        if res["tec_data"][sample][target]["Ncopy_mean"] > 0.0:
                            if res["ref_data"][sample]["Ncopy_gem"] > 0.0:
                                relEx = res["tec_data"][sample][target]["Ncopy_mean"] / (res["ref_data"][sample]["Ncopy_gem"] / res["ref_geoAll"])
                        res["rel_data"][sample][target]["rel_expression"] = relEx

            for sample in res["rel_data"]:
                for target in res["rel_data"][sample]:
                    res["rel_data"][sample][target]["raw_vals"] = []
                    if res["rel_data"][sample][target]["rel_expression"] > 0.0:
                        for indivVal in nCopyData["Ncopy"][sample][target]:
                            res["rel_data"][sample][target]["raw_vals"].append((indivVal / res["ref_data"][sample]["Ncopy_gem"]))

        if overlapType == "annotation":
            res["anno_data"] = {}
            res["anno_stats"] = {}
            res["anno_key"] = selAnnotation
            for sample in res["rel_data"]:
                for target in res["rel_data"][sample]:
                    if len(selReferences) > 0:
                        if res["rel_data"][sample][target]["rel_expression"] > 0.0:
                            if sample in samSelAnno:
                                if target not in res["anno_data"]:
                                    res["anno_data"][target] = {}
                                if samSelAnno[sample] not in res["anno_data"][target]:
                                    res["anno_data"][target][samSelAnno[sample]] = {}
                                if "raw_vals" not in res["anno_data"][target][samSelAnno[sample]]:
                                    res["anno_data"][target][samSelAnno[sample]]["raw_vals"] = []
                                res["anno_data"][target][samSelAnno[sample]]["raw_vals"].append(res["rel_data"][sample][target]["rel_expression"])
                    else:
                        if res["tec_data"][sample][target]["Ncopy_mean"] > 0.0:
                            if sample in samSelAnno:
                                if target not in res["anno_data"]:
                                    res["anno_data"][target] = {}
                                if samSelAnno[sample] not in res["anno_data"][target]:
                                    res["anno_data"][target][samSelAnno[sample]] = {}
                                if "raw_vals" not in res["anno_data"][target][samSelAnno[sample]]:
                                    res["anno_data"][target][samSelAnno[sample]]["raw_vals"] = []
                                res["anno_data"][target][samSelAnno[sample]]["raw_vals"].append(res["tec_data"][sample][target]["Ncopy_mean"])
            for target in res["anno_data"]:
                statTarGroup = []
                translateGrp = {}
                tarCountUp = 0
                for annoVal in res["anno_data"][target]:
                    annoCollVals = res["anno_data"][target][annoVal]["raw_vals"]
                    statTarGroup.append(annoCollVals)
                    translateGrp[tarCountUp] = annoVal
                    tarCountUp += 1
                    if len(annoCollVals) == 0:
                        res["anno_data"][target][annoVal]["mean"] = -1.0
                        res["anno_data"][target][annoVal]["sem"] = -1.0
                    else:
                        res["anno_data"][target][annoVal]["mean"] = float(np.mean(annoCollVals))
                        if len(annoCollVals) == 1:
                            res["anno_data"][target][annoVal]["sem"] = 0.0
                        else:
                            res["anno_data"][target][annoVal]["sem"] = float(scp.sem(annoCollVals))
                res["anno_stats"][target] = runStatistics(statTarGroup, statsParametric, translateGrp, statAlpha)

        if saveResultsCSV:
            if len(res["dil_std"]) > 0:
                res["tsv"]["dil_std"] = "Target\tCurve PCR Eff.\tDilution PCR Eff.\tStd Max\tNcopy Max\tSTD Min\tNcopy Min\tStd Ratio\tNcopy Ratio\t"
                res["tsv"]["dil_std"] += "Slope Bias\tCorrelation R\tLinearity\tPrecision\tResolution\n"
                sortTars = sorted(res["dil_std"].keys())
                for target in sortTars:
                    res["tsv"]["dil_std"] += target + "\t"
                    res["tsv"]["dil_std"] += "{:.4f}".format(res["dil_std"][target]["curve_pcr_eff"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.4f}".format(res["dil_std"][target]["dilution_pcr_eff"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["expected_max"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["mean_max"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["expected_min"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["mean_min"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["expected_ratio"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.2f}".format(res["dil_std"][target]["mean_ratio"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.6f}".format(res["dil_std"][target]["slope_bias"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.6f}".format(res["dil_std"][target]["correlation_R"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.6f}".format(res["dil_std"][target]["linearity"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.6f}".format(res["dil_std"][target]["reproducibility"]) + "\t"
                    res["tsv"]["dil_std"] += "{:.6f}".format(res["dil_std"][target]["detectable_diff"]) + "\t"
                    res["tsv"]["dil_std"] += "\n"

            res["tsv"]["technical_data"] = "Sample\tSample Type\t"
            if inclAnnotation:
                for currAnno in sortedAnnoKeys:
                    res["tsv"]["technical_data"] += currAnno + "\t"
            res["tsv"]["technical_data"] += "Target\tTarget Type\tError\tNote\t"
            res["tsv"]["technical_data"] += "n Tec. Rep.\tMean Ncopy\tSD Ncopy\tCV Ncopy\tIndividual Values\n"
            sortSam = sorted(res["tec_data"].keys())
            for sample in sortSam:
                sortTar = sorted(res["tec_data"][sample].keys())
                for target in sortTar:
                    res["tsv"]["technical_data"] += sample + "\t"
                    res["tsv"]["technical_data"] += res["tec_data"][sample][target]["sample_type"] + "\t"
                    if inclAnnotation:
                        for currAnno in sortedAnnoKeys:
                            annoVal = ""
                            if sample in samAllAnnos:
                                if currAnno in samAllAnnos[sample]:
                                    annoVal = samAllAnnos[sample][currAnno]
                            res["tsv"]["technical_data"] += annoVal + "\t"
                    res["tsv"]["technical_data"] += target + "\t"
                    res["tsv"]["technical_data"] += res["tec_data"][sample][target]["target_type"] + "\t"
                    res["tsv"]["technical_data"] += res["tec_data"][sample][target]["error"] + "\t"
                    res["tsv"]["technical_data"] += res["tec_data"][sample][target]["note"] + "\t"
                    res["tsv"]["technical_data"] += str(res["tec_data"][sample][target]["n_tec_rep"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.4f}".format(res["tec_data"][sample][target]["Ncopy_mean"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.4f}".format(res["tec_data"][sample][target]["Ncopy_sd"]) + "\t"
                    res["tsv"]["technical_data"] += "{:.6f}".format(res["tec_data"][sample][target]["Ncopy_cv"]) + "\t"
                    for indivVal in res["tec_data"][sample][target]["raw_vals"]:
                        res["tsv"]["technical_data"] += "{:.4f}".format(indivVal) + ";"
                    res["tsv"]["technical_data"] += "\n"

            if len(selReferences) > 0:
                res["tsv"]["reference_data"] = "Sample\t"
                if inclAnnotation:
                    for currAnno in sortedAnnoKeys:
                        res["tsv"]["reference_data"] += currAnno + "\t"
                res["tsv"]["reference_data"] += "Error\tn Ref. Genes\tGeometric Mean Ncopy\tRelative Expresssion\n"
                sortSam = sorted(res["ref_data"].keys())
                for sample in sortSam:
                    res["tsv"]["reference_data"] += sample + "\t"
                    if inclAnnotation:
                        for currAnno in sortedAnnoKeys:
                            annoVal = ""
                            if sample in samAllAnnos:
                                if currAnno in samAllAnnos[sample]:
                                    annoVal = samAllAnnos[sample][currAnno]
                            res["tsv"]["reference_data"] += annoVal + "\t"
                    if res["ref_data"][sample]["ref_missing"]:
                        res["tsv"]["reference_data"] += "Reference Genes without Ncopy"
                    res["tsv"]["reference_data"] += "\t"
                    res["tsv"]["reference_data"] += str(res["ref_data"][sample]["Ncopy_num"]) + "\t"
                    res["tsv"]["reference_data"] += "{:.4f}".format(res["ref_data"][sample]["Ncopy_gem"]) + "\t"
                    res["tsv"]["reference_data"] += "{:.4f}".format(res["ref_data"][sample]["Ncopy_gem"] / res["ref_geoAll"]) + "\t"
                    # for indivVal in res["ref_data"][sample]["raw_vals"]:
                    #     res["tsv"]["reference_data"] += "{:.4e}".format(indivVal) + ";"
                    res["tsv"]["reference_data"] += "\n"
                res["tsv"]["reference_data"] += "Geometric Mean\t\t\t\t" + "{:.4f}".format(res["ref_geoAll"]) + "\t1.0000\n"

                res["tsv"]["relative_data"] = "Sample\tSample Type\t"
                if inclAnnotation:
                    for currAnno in sortedAnnoKeys:
                        res["tsv"]["relative_data"] += currAnno + "\t"
                res["tsv"]["relative_data"] += "Target\tTarget Type\tn Rel. Exp.\tRel. Expression\n"
                sortSam = sorted(res["rel_data"].keys())
                for sample in sortSam:
                    sortTar = sorted(res["rel_data"][sample].keys())
                    for target in sortTar:
                        res["tsv"]["relative_data"] += sample + "\t"
                        res["tsv"]["relative_data"] += res["tec_data"][sample][target]["sample_type"] + "\t"
                        if inclAnnotation:
                            for currAnno in sortedAnnoKeys:
                                annoVal = ""
                                if sample in samAllAnnos:
                                    if currAnno in samAllAnnos[sample]:
                                        annoVal = samAllAnnos[sample][currAnno]
                                res["tsv"]["relative_data"] += annoVal + "\t"
                        res["tsv"]["relative_data"] += target + "\t"
                        res["tsv"]["relative_data"] += res["tec_data"][sample][target]["target_type"] + "\t"
                        res["tsv"]["relative_data"] += str(len(res["rel_data"][sample][target]["raw_vals"])) + "\t"
                        if res["rel_data"][sample][target]["ref_missing"]:
                            res["tsv"]["relative_data"] += "Reference Genes without Ncopy"
                        res["tsv"]["relative_data"] += "{:.4f}".format(res["rel_data"][sample][target]["rel_expression"]) + "\t"
                        # for indivVal in res["rel_data"][sample][target]["raw_vals"]:
                        #     res["tsv"]["relative_data"] += "{:.6f}".format(indivVal) + ";"
                        res["tsv"]["relative_data"] += "\n"

            if overlapType == "annotation":
                res["tsv"]["annotation_data"] = res["anno_key"] + "\t"
                if inclAnnotation:
                    for currAnno in sortedAnnoKeys:
                        if currAnno == selAnnotation:
                            continue
                        if currAnno not in overSelAnno:
                            continue
                        if overSelAnno[currAnno]["conf"]:
                            continue
                        res["tsv"]["annotation_data"] += currAnno + "\t"
                res["tsv"]["annotation_data"] += "Target\tn\tRel. Expression\tSEM\tIndividual Values\n"
                sortTar = sorted(res["anno_data"].keys())
                for target in sortTar:
                    sortAnno = sorted(res["anno_data"][target].keys())
                    for annoVal in sortAnno:
                        res["tsv"]["annotation_data"] += annoVal + "\t"
                        if inclAnnotation:
                            for currAnno in sortedAnnoKeys:
                                if currAnno == selAnnotation:
                                    continue
                                if currAnno not in overSelAnno:
                                    continue
                                if overSelAnno[currAnno]["conf"]:
                                    continue
                                annoData = ""
                                if currAnno in overSelAnno:
                                    if annoVal in overSelAnno[currAnno]["data"]:
                                        annoData = overSelAnno[currAnno]["data"][annoVal]
                                res["tsv"]["annotation_data"] += annoData + "\t"
                        res["tsv"]["annotation_data"] += target + "\t"
                        res["tsv"]["annotation_data"] += str(len(res["anno_data"][target][annoVal]["raw_vals"])) + "\t"
                        res["tsv"]["annotation_data"] += "{:.4f}".format(res["anno_data"][target][annoVal]["mean"]) + "\t"
                        res["tsv"]["annotation_data"] += "{:.4f}".format(res["anno_data"][target][annoVal]["sem"]) + "\t"
                        for indivVal in res["anno_data"][target][annoVal]["raw_vals"]:
                            res["tsv"]["annotation_data"] += "{:.4f}".format(indivVal) + ";"
                        res["tsv"]["annotation_data"] += "\n"

                res["tsv"]["statistics_data"] = "Target\tStat. Test\tAlpha\tTest Var. Name\tTest Output\tdf\tP Value\n"
                for target in sortTar:
                    res["tsv"]["statistics_data"] += target + "\t"
                    res["tsv"]["statistics_data"] += res["anno_stats"][target]["test name"] + "\t"
                    res["tsv"]["statistics_data"] += str(res["anno_stats"][target]["alpha"]) + "\t"
                    res["tsv"]["statistics_data"] += res["anno_stats"][target]["stat name"] + "\t"
                    res["tsv"]["statistics_data"] += "{:.4f}".format(res["anno_stats"][target]["stat val"]) + "\t"
                    res["tsv"]["statistics_data"] += str(res["anno_stats"][target]["df"]) + "\t"
                    res["tsv"]["statistics_data"] += "{:.4f}".format(res["anno_stats"][target]["p val"]) + "\n"

                res["tsv"]["statistics_multi_comp"] = ""
                for target in sortTar:
                    if res["anno_stats"][target]["multi comparison"]  != "":
                        res["tsv"]["statistics_multi_comp"] += target + "\n"
                        for row in range(0, len(res["anno_stats"][target]["multi comparison"])):
                            for col in range(0, len(res["anno_stats"][target]["multi comparison"][row])):
                                res["tsv"]["statistics_multi_comp"] += str(res["anno_stats"][target]["multi comparison"][row][col]) + "\t"
                            res["tsv"]["statistics_multi_comp"] = re.sub(r"\t$", "\n", res["tsv"]["statistics_multi_comp"])
                        res["tsv"]["statistics_multi_comp"] += "\n\n"

        return res


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

    def remove_calculations(self):
        """Removes all values which result of calculations in a run.

        Args:
            self: The class self parameter.

        Returns:
            A string with the data.
        """

        ret = []
        hint1 = ""
        hint2 = ""
        hint3 = ""
        hint4 = ""
        hint5 = ""
        hint6 = ""
        hint7 = ""
        hint8 = ""
        hint9 = ""
        hint10 = ""
        hint11 = ""
        hint12 = ""
        hint13 = ""
        hint14 = ""
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                exp5 = _get_all_children(react_data, "cq")
                for node5 in exp5:
                    hint1 = "Deleted react data \"cq\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "N0")
                for node5 in exp5:
                    hint2 = "Deleted react data \"N0\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "Ncopy")
                for node5 in exp5:
                    hint14 = "Deleted react data \"Ncopy\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "ampEffMet")
                for node5 in exp5:
                    hint3 = "Deleted react data \"ampEffMet\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "ampEff")
                for node5 in exp5:
                    hint4 = "Deleted react data \"ampEff\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "ampEffSE")
                for node5 in exp5:
                    hint5 = "Deleted react data \"ampEffSE\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "corrF")
                for node5 in exp5:
                    hint6 = "Deleted react data \"corrF\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "corrP")
                for node5 in exp5:
                    hint7 = "Deleted react data \"corrP\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "corrCq")
                for node5 in exp5:
                    hint8 = "Deleted react data \"corrCq\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "meltTemp")
                for node5 in exp5:
                    hint9 = "Deleted react data \"meltTemp\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "endPt")
                for node5 in exp5:
                    hint10 = "Deleted react data \"endPt\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "bgFluor")
                for node5 in exp5:
                    hint11 = "Deleted react data \"bgFluor\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "bgFluorSlp")
                for node5 in exp5:
                    hint12 = "Deleted react data \"bgFluorSlp\" elements."
                    react_data.remove(node5)
                exp5 = _get_all_children(react_data, "quantFluor")
                for node5 in exp5:
                    hint13 = "Deleted react data \"quantFluor\" elements."
                    react_data.remove(node5)

        if hint1 != "":
            ret.append(hint1)
        if hint2 != "":
            ret.append(hint2)
        if hint14 != "":
            ret.append(hint14)
        if hint3 != "":
            ret.append(hint3)
        if hint4 != "":
            ret.append(hint4)
        if hint5 != "":
            ret.append(hint5)
        if hint6 != "":
            ret.append(hint6)
        if hint7 != "":
            ret.append(hint7)
        if hint8 != "":
            ret.append(hint8)
        if hint9 != "":
            ret.append(hint9)
        if hint10 != "":
            ret.append(hint10)
        if hint11 != "":
            ret.append(hint11)
        if hint12 != "":
            ret.append(hint12)
        if hint13 != "":
            ret.append(hint13)
        return ret

    def remove_notes(self):
        """Removes all values which result of calculations in a run.

        Args:
            self: The class self parameter.

        Returns:
            A string with the data.
        """

        ret = []
        hint1 = ""
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                exp5 = _get_all_children(react_data, "note")
                for node5 in exp5:
                    hint1 = "Deleted react data \"note\" elements."
                    react_data.remove(node5)

        if hint1 != "":
            ret.append(hint1)
        return ret

    def remove_errors(self):
        """Removes all values which result of calculations in a run.

        Args:
            self: The class self parameter.

        Returns:
            A string with the data.
        """

        ret = []
        hint1 = ""
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                exp5 = _get_all_children(react_data, "excl")
                for node5 in exp5:
                    hint1 = "Deleted react data \"error\" elements."
                    react_data.remove(node5)

        if hint1 != "":
            ret.append(hint1)
        return ret

    def export_table(self, dMode):
        """Returns a tab seperated table file with the react fluorescence data in RDES format.

        Args:
            self: The class self parameter.
            dMode: amp for amplification data, melt for meltcurve data

        Returns:
            A string with the data.
        """

        tarTypeLookup = {}
        tarDyeLookup = {}
        data = ""

        # Get the information for the lookup dictionaries
        pExp = self._node.getparent()
        pRoot = pExp.getparent()
        transSamTar = _sampleTypeToDics(pRoot)
        targets = _get_all_children(pRoot, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType != "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None:
                    if forId.attrib['id'] != "":
                        tarDyeLookup[tarId] = forId.attrib['id']

        # Now create the header line
        data += "Well\tSample\tSample Type\tTarget\tTarget Type\tDye\t"
        if dMode == "amp":
            data += "Cq"
        else:
            data += "Tm"
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
            headArr = sorted(headArr, key=float)
        else:
            mdps = _get_all_children(react_datas[0], "mdp")
            for mdp in mdps:
                headArr.append(_get_first_child_text(mdp, "tmp"))
            headArr = sorted(headArr, key=float)
        for hElem in headArr:
            data += "\t" + hElem
        data += '\n'

        # Now create the data lines
        reacts = _get_all_children(self._node, "react")
        wellData = []
        for react in reacts:
            reactId = react.get('id')
            pWell = str(reactId)
            if int(self["pcrFormat_columns"]) != 1 and int(self["pcrFormat_rows"]) != 1:
                pIdNumber = (int(reactId) - 1) % int(self["pcrFormat_columns"]) + 1
                pIdLetter = chr(ord("A") + int((int(reactId) - 1) / int(self["pcrFormat_columns"])))
                pWell = pIdLetter + str(pIdNumber)
            dataSample = pWell + '\t'
            react_sample = "No Sample"
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    react_sample = forId.attrib['id']
            dataSample += react_sample + '\t'
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
                react_sample_type = "No Sample Type"
                if react_sample in transSamTar:
                    if react_target in transSamTar[react_sample]:
                        react_sample_type = transSamTar[react_sample][react_target]
                dataLine = dataSample + react_sample_type
                dataLine += "\t" + react_target + '\t' + react_target_type + '\t' + react_target_dye
                col7 = ""
                if dMode == "amp":
                    col7 = _get_first_child_text(react_data, "cq")
                else:
                    col7 = _get_first_child_text(react_data, "meltTemp")
                dataLine += '\t' + col7
                fluorList = []
                if dMode == "amp":
                    adps = _get_all_children(react_data, "adp")
                    for adp in adps:
                        cyc = _get_first_child_text(adp, "cyc")
                        fluor = _get_first_child_text(adp, "fluor")
                        fluorList.append([cyc, fluor])
                    fluorList = sorted(fluorList, key=_sort_list_float)
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
            dMode: amp for amplification data, melt for meltcurve data or ncopy for only Ncopy data.

        Returns:
            A string with the modifications made.
        """

        dataVersion = rootEl._node.attrib['version'] # _get_first_child_text(rootEl, "version")
        ret = ""
        negWarning = False
        with open(filename, "r") as tfile:
            fileContent = tfile.read()

        newlineFix = fileContent.replace("\r\n", "\n")
        tabLines = newlineFix.split("\n")

        head = tabLines[0].split("\t")
        if (head[0].lower() != "well" or head[1].lower() != "sample" or head[2].lower() != "sample type" or
                head[3].lower() != "target" or head[4].lower() != "target type" or head[5].lower() != "dye"):
            raise RdmlError('The tab-format is not valid, essential columns are missing.')

        if head[6].lower() in ["cq", "ncopy"]:
            if dMode == "melt":
                raise RdmlError('Import Format Error: Amplification data immported as melting data.')
        if head[6].lower() == "tm":
            if dMode != "melt":
                raise RdmlError('Import Format Error: Melting data immported as amplification data.')

        # Get the information for the lookup dictionaries
        samTypeLookup = {}
        tarTypeLookup = {}
        dyeLookup = {}
        transSamTar = _sampleTypeToDics(rootEl._node)
        samples = _get_all_children(rootEl._node, "sample")
        for sample in samples:
            if sample.attrib['id'] != "":
                samId = sample.attrib['id']
                samTypeLookup[samId] = "unkn"
                forType = _get_all_children(sample, "type")
                for subNode in forType:
                    if 'targetId' not in subNode.attrib:
                        samTypeLookup[samId] = subNode.text
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType != "":
                    tarTypeLookup[tarId] = forType
                forId = _get_first_child(target, "dyeId")
                if forId is not None and forId.attrib['id'] != "":
                    dyeLookup[forId.attrib['id']] = 1

        # Process the lines
        for tabLine in tabLines[1:]:
            sLin = tabLine.split("\t")
            if (len(sLin) < 7 or sLin[0] == "" or sLin[1] == ""
                              or sLin[3] == "" or sLin[5] == ""):
                ret += "Skipped reaction \"" + sLin[0] + "\"\n"
                continue

            # Undefined are set to "unkn" and "toi"
            if sLin[2] not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
                sLin[2] = "unkn"
            if sLin[4] not in ["ref", "toi"]:
                sLin[4] = "toi"

            if sLin[1] not in samTypeLookup:
                rootEl.new_sample(sLin[1])
                samEl = rootEl.get_sample(byid=sLin[1])
                samEl.new_type(sLin[2])
                samTypeLookup[sLin[1]] = sLin[2]
                ret += "Created sample \"" + sLin[1] + "\" with type \"" + sLin[2] + "\"\n"
                if sLin[1] not in transSamTar:
                    transSamTar[sLin[1]] = {}
            if sLin[2] != samTypeLookup[sLin[1]]:
                if ((sLin[3] not in transSamTar[sLin[1]]) or
                        (transSamTar[sLin[1]][sLin[3]] != sLin[2])):
                    samEl = rootEl.get_sample(byid=sLin[1])
                    samEl.new_type(sLin[2], targetId=sLin[3])
                    transSamTar[sLin[1]][sLin[3]] = sLin[2]
                    ret += "Created sample \"" + sLin[1] + "\" with type \"" + sLin[2] + "\" for target \"" + sLin[3] + "\"\n"
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
            if re.search(r"\D\d+", sLin[0]):
                old_letter = ord(re.sub(r"\d", "", sLin[0]).upper()) - ord("A")
                old_nr = int(re.sub(r"\D", "", sLin[0]))
                newId = old_nr + old_letter * int(self["pcrFormat_columns"])
                wellPos = str(newId)
            if re.search(r"\D\d+\D\d+", sLin[0]):
                old_left = re.sub(r"\D\d+$", "", sLin[0])
                old_left_letter = ord(re.sub(r"\d", "", old_left).upper()) - ord("A")
                old_left_nr = int(re.sub(r"\D", "", old_left)) - 1
                old_right = re.sub(r"^\D\d+", "", sLin[0])
                old_right_letter = ord(re.sub(r"\d", "", old_right).upper()) - ord("A")
                old_right_nr = int(re.sub(r"\D", "", old_right))
                newId = old_left_nr * 8 + old_right_nr + old_left_letter * 768 + old_right_letter * 96
                wellPos = str(newId)

            exp = _get_all_children(self._node, "react")
            for node in exp:
                if wellPos == node.attrib['id']:
                    react = node
                    forId = _get_first_child_text(react, "sample")
                    if forId and forId != "" and forId.attrib['id'] != sLin[1]:
                        ret += "Missmatch: Well " + wellPos + " (" + sLin[0] + ") has sample \"" + forId.attrib['id'] + \
                               "\" in RDML file and sample \"" + sLin[1] + "\" in tab file.\n"
                    break
            if react is None:
                new_node = et.Element("react", id=wellPos)
                place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                self._node.insert(place, new_node)
                react = new_node
                new_node = et.Element("sample", id=sLin[1])
                react.insert(0, new_node)

            exp = _get_all_children(react, "data")
            for node in exp:
                forId = _get_first_child(node, "tar")
                if forId is not None and forId.attrib['id'] == sLin[3]:
                    data = node
                    break
            if data is None:
                new_node = et.Element("data")
                place = _get_tag_pos(react, "data", ["sample", "data", "partitions"], 9999999)
                react.insert(place, new_node)
                data = new_node
                new_node = et.Element("tar", id=sLin[3])
                place = _get_tag_pos(data, "tar",
                                     _getXMLDataType(),
                                     9999999)
                data.insert(place, new_node)

            keyFor7 = "cq"
            if dMode == "melt":
                keyFor7 = "meltTemp"
            if head[6].lower() == "ncopy":
                keyFor7 = "Ncopy"

            present7 = _get_first_child(data, keyFor7)
            cont7 = re.sub(r'[^eE0-9\.;\-]', '', sLin[6])
            if cont7 == "" or (cont7 is None):
                if present7 is not None:
                    data.remove(present7)
            else:
                if ((dMode == "melt" and dataVersion in ["1.3", "1.4"]) or 
                    (dMode in ["amp", "ncopy"] and dataVersion in ["1.4"] and keyFor7 == "Ncopy") or
                    (dMode == "amp" and keyFor7 == "cq")):
                    if present7 is None:
                        new_node = et.Element(keyFor7)
                        new_node.text = cont7
                        place = _get_tag_pos(data, keyFor7,
                                            _getXMLDataType(),
                                            9999999)
                        data.insert(place, new_node)
                    else:
                        present7.text = cont7

            if dMode == "amp":
                presentAmp = _get_first_child(data, "adp")
                if presentAmp is not None:
                    ret += "Well " + wellPos + " (" + sLin[0] + ") with sample \"" + sLin[1] + " and target \"" + \
                           sLin[3] + "\" has already amplification data, no data were added.\n"
                else:
                    colCount = 7
                    uniqueCheck = []
                    for col in sLin[7:]:
                        if head[colCount] not in uniqueCheck:
                            uniqueCheck.append(head[colCount])
                            try:
                                colToFloat = float(col)
                            except ValueError:
                                pass
                            else:
                                if math.isfinite(colToFloat):
                                    if colToFloat < 0.0:
                                        negWarning = True
                                    new_node = et.Element("adp")
                                    place = _get_tag_pos(data, "adp",
                                                        _getXMLDataType(),
                                                        9999999)
                                    data.insert(place, new_node)
                                    new_sub = et.Element("cyc")
                                    new_sub.text = head[colCount]
                                    place = _get_tag_pos(new_node, "cyc", ["cyc", "tmp", "fluor"], 9999999)
                                    new_node.insert(place, new_sub)
                                    new_sub = et.Element("fluor")
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
                    colCount = 7
                    uniqueCheck = []
                    for col in sLin[7:]:
                        if head[colCount] not in uniqueCheck:
                            uniqueCheck.append(head[colCount])
                            try:
                                colToFloat = float(col)
                            except ValueError:
                                pass
                            else:
                                if math.isfinite(colToFloat):
                                    if colToFloat < 0.0:
                                        negWarning = True
                                    new_node = et.Element("mdp")
                                    place = _get_tag_pos(data, "mdp",
                                                        _getXMLDataType(),
                                                        9999999)
                                    data.insert(place, new_node)
                                    new_sub = et.Element("tmp")
                                    new_sub.text = head[colCount]
                                    place = _get_tag_pos(new_node, "tmp", ["tmp", "fluor"], 9999999)
                                    new_node.insert(place, new_sub)
                                    new_sub = et.Element("fluor")
                                    new_sub.text = col
                                    place = _get_tag_pos(new_node, "fluor", ["tmp", "fluor"], 9999999)
                                    new_node.insert(place, new_sub)
                        colCount += 1
        if negWarning == True:
            if ret == "":
                ret = "Warning: Negative fluorescence values detected. Be sure to import raw data without baseline correction.\n"
            else:
                ret = "Warning: Negative fluorescence values detected. Be sure to import raw data without baseline correction.\n" + ret
                ret += "\nWarning: Negative fluorescence values detected. Be sure to import raw data without baseline correction."
        return ret

    def import_digital_data(self, rootEl, fileformat, filename, filelist, ignoreCh=""):
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

        partElemLS = _getXMLPartitionDataType()
        tempList = re.split(r"\D+", ignoreCh)
        ignoreList = []
        for posNum in tempList:
            if re.search(r"\d", posNum):
                ignoreList.append(int(posNum))

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
                samTypeLookup[samId] = "unkn"
                forType = _get_first_child_text(sample, "type")
                if forType != "":
                    samTypeLookup[samId] = forType
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType != "":
                    tarTypeLookup[tarId] = forType
        dyes = _get_all_children(rootEl._node, "dye")
        for dye in dyes:
            if dye.attrib['id'] != "":
                dyeLookup[dye.attrib['id']] = 1
        # Work the overview file
        if filename is not None:
            with open(filename, newline='') as tfile:  # add encoding='utf-8' ?
                posCount = 0
                posWell = 0
                posSample = -1
                posSampleType = -1
                posDye = -1
                posDyeCh2 = -1
                posDyeCh3 = -1
                posTarget = -1
                posTargetCh2 = -1
                posTargetCh3 = -1
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

                countUpTarget = 1

                if fileformat == "RDML":
                    tabLines = list(csv.reader(tfile, delimiter='\t'))
                    for hInfo in tabLines[0]:
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
                    tabLines = list(csv.reader(tfile, delimiter=','))
                    for hInfo in tabLines[0]:
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
                    tabLines = list(csv.reader(tfile, delimiter=','))
                    for hInfo in tabLines[0]:
                        hInfo = re.sub(r"^ +", '', hInfo)
                        if hInfo == "SampleName":
                            posSample = posCount
                        # This is a hack of the format to allow specification of targets
                        if hInfo == "Blue_Channel_Target":
                            posTarget = posCount
                        if hInfo == "Green_Channel_Target":
                            posTargetCh2 = posCount
                        if hInfo == "Red_Channel_Target":
                            posTargetCh3 = posCount
                        # End of hack
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
                for rowNr in range(1, len(tabLines)):
                    emptyLine = True
                    if len(tabLines[rowNr]) < 7:
                        continue
                    for colNr in range(0, len(tabLines[rowNr])):
                        if tabLines[rowNr][colNr] != "":
                            emptyLine = False
                            tabLines[rowNr][colNr] = re.sub(r'^ +', '', tabLines[rowNr][colNr])
                            tabLines[rowNr][colNr] = re.sub(r' +$', '', tabLines[rowNr][colNr])
                    if emptyLine is True:
                        continue
                    sLin = tabLines[rowNr]

                    if sLin[posSample] not in samTypeLookup:
                        posSampleTypeName = "unkn"
                        if posSampleType != -1:
                            posSampleTypeName = sLin[posSampleType]
                        rootEl.new_sample(sLin[posSample])
                        samEl = rootEl.get_sample(byid=sLin[posSample])
                        samEl.new_type(sLin[posSampleTypeName])
                        samTypeLookup[sLin[posSample]] = posSampleTypeName
                        ret += "Created sample \"" + sLin[posSample] + "\" with type \"" + posSampleTypeName + "\"\n"

                    # Fix well position
                    wellPos = re.sub(r"\"", "", sLin[posWell])
                    if fileformat == "Stilla":
                        wellPos = re.sub(r'^\d+-', '', wellPos)

                    # Create nonexisting targets and dyes
                    if fileformat == "Stilla":
                        if 1 not in ignoreList:
                            if posTarget > -1:
                                crTarName = sLin[posTarget]
                            else:
                                crTarName = " Target " + str(countUpTarget) + " Ch1"
                                countUpTarget += 1
                            chan = "Ch1"
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
                            if wellPos.upper() not in headerLookup:
                                headerLookup[wellPos.upper()] = {}
                            headerLookup[wellPos.upper()][chan] = crTarName
                        if 2 not in ignoreList:
                            if posTargetCh2 > -1:
                                crTarName = sLin[posTargetCh2]
                            else:
                                crTarName = " Target " + str(countUpTarget) + " Ch2"
                                countUpTarget += 1
                            chan = "Ch2"
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
                            if wellPos.upper() not in headerLookup:
                                headerLookup[wellPos.upper()] = {}
                            headerLookup[wellPos.upper()][chan] = crTarName
                        if 3 not in ignoreList:
                            if posTargetCh3 > -1:
                                crTarName = sLin[posTargetCh3]
                            else:
                                crTarName = " Target " + str(countUpTarget) + " Ch3"
                                countUpTarget += 1
                            chan = "Ch3"
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
                            if wellPos.upper() not in headerLookup:
                                headerLookup[wellPos.upper()] = {}
                            headerLookup[wellPos.upper()][chan] = crTarName
                    else:
                        if fileformat == "Bio-Rad":
                            posDyeName = sLin[posDye][:3]
                        else:
                            posDyeName = sLin[posDye]
                        if posTarget > -1 and int(re.sub(r"\D", "", posDyeName)) not in ignoreList:
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

                            if wellPos.upper() not in headerLookup:
                                headerLookup[wellPos.upper()] = {}
                            headerLookup[wellPos.upper()][posDyeName] = sLin[posTarget]

                    if posFilename != -1 and sLin[posFilename] != "":
                        fileNameSuggLookup[wellPos.upper()] = sLin[posFilename]

                    react = None
                    partit = None
                    data = None

                    # Get the position number if required
                    wellPosStore = wellPos
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
                            if forId and forId != "" and forId.attrib['id'] != sLin[posSample]:
                                ret += "Missmatch: Well " + wellPos + " (" + sLin[posWell] + ") has sample \"" + forId.attrib['id'] + \
                                       "\" in RDML file and sample \"" + sLin[posSample] + "\" in tab file.\n"
                            break
                    if react is None:
                        new_node = et.Element("react", id=wellPos)
                        place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                        self._node.insert(place, new_node)
                        react = new_node
                        new_node = et.Element("sample", id=sLin[posSample])
                        react.insert(0, new_node)

                    partit = _get_first_child(react, "partitions")
                    if partit is None:
                        new_node = et.Element("partitions")
                        place = _get_tag_pos(react, "partitions", ["sample", "data", "partitions"], 9999999)
                        react.insert(place, new_node)
                        partit = new_node
                        new_node = et.Element("volume")
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
                            if i in ignoreList:
                                continue
                            data = None
                            posDyeName = "Ch" + str(i)
                            stillaTarget = headerLookup[wellPosStore.upper()][posDyeName]
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

                            for node in exp:
                                forId = _get_first_child(node, "tar")
                                if forId is not None and forId.attrib['id'] == stillaTarget:
                                    data = node
                                    break

                            if data is None:
                                new_node = et.Element("data")
                                place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                                partit.insert(place, new_node)
                                data = new_node
                                new_node = et.Element("tar", id=stillaTarget)
                                place = _get_tag_pos(data, "tar", partElemLS, 9999999)
                                data.insert(place, new_node)

                            new_node = et.Element("pos")
                            new_node.text = stillaPos
                            place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                            data.insert(place, new_node)

                            new_node = et.Element("neg")
                            new_node.text = stillaNeg
                            place = _get_tag_pos(data, "neg", partElemLS, 9999999)
                            data.insert(place, new_node)

                            new_node = et.Element("conc")
                            new_node.text = stillaConc
                            place = _get_tag_pos(data, "conc", partElemLS, 9999999)
                            data.insert(place, new_node)
                    else:
                        exp = _get_all_children(partit, "data")
                        for node in exp:
                            forId = _get_first_child(node, "tar")
                            if forId is not None and forId.attrib['id'] == sLin[posTarget]:
                                data = node
                                break
                        if data is None:
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            data = new_node
                            new_node = et.Element("tar", id=sLin[posTarget])
                            place = _get_tag_pos(data, "tar", partElemLS, 9999999)
                            data.insert(place, new_node)

                        new_node = et.Element("pos")
                        new_node.text = sLin[posPositives]
                        place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                        data.insert(place, new_node)

                        new_node = et.Element("neg")
                        new_node.text = sLin[posNegatives]
                        place = _get_tag_pos(data, "neg", partElemLS, 9999999)
                        data.insert(place, new_node)

                        if posUndefined != -1 and sLin[posUndefined] != "":
                            new_node = et.Element("undef")
                            new_node.text = sLin[posUndefined]
                            place = _get_tag_pos(data, "neg", partElemLS, 9999999)
                            data.insert(place, new_node)

                        if posExcluded != -1 and sLin[posExcluded] != "":
                            new_node = et.Element("excl")
                            new_node.text = sLin[posExcluded]
                            place = _get_tag_pos(data, "neg", partElemLS, 9999999)
                            data.insert(place, new_node)

                        if posCopConc != -1:
                            new_node = et.Element("conc")
                            if int(sLin[posPositives]) == 0:
                                new_node.text = "0"
                            else:
                                if fileformat == "RDML":
                                    new_node.text = sLin[posCopConc]
                                elif fileformat == "Bio-Rad":
                                    new_node.text = str(float(sLin[posCopConc])/20)
                                else:
                                    new_node.text = sLin[posCopConc]
                            place = _get_tag_pos(data, "conc", partElemLS, 9999999)
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
            if re.search(r"\D\d+", well):
                old_letter = ord(re.sub(r"\d", "", well).upper()) - ord("A")
                old_nr = int(re.sub(r"\D", "", well))
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
                    rootEl.new_sample(sampleName)
                    samEl = rootEl.get_sample(byid=sampleName)
                    samEl.new_type("unkn")
                    samTypeLookup[sampleName] = "unkn"
                    ret += "Created sample \"" + sampleName + "\" with type \"" + "unkn" + "\"\n"
                new_node = et.Element("react", id=wellPos)
                place = _get_tag_pos(self._node, "react", self.xmlkeys(), 9999999)
                self._node.insert(place, new_node)
                react = new_node
                new_node = et.Element("sample", id=sampleName)
                react.insert(0, new_node)

            partit = _get_first_child(react, "partitions")
            if partit is None:
                new_node = et.Element("partitions")
                place = _get_tag_pos(react, "partitions", ["sample", "data", "partitions"], 9999999)
                react.insert(place, new_node)
                partit = new_node
                new_node = et.Element("volume")
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

            # print(finalFileName, flush=True)

            with open(fileLookup[well], newline='') as wellfile:  # add encoding='utf-8' ?
                if fileformat == "RDML":
                    wellLines = list(csv.reader(wellfile, delimiter='\t'))
                    wellFileContent = wellfile.read()
                    _writeFileInRDML(self._rdmlFilename, finalFileName, wellFileContent)

                    delElem = _get_first_child(partit, "endPtTable")
                    if delElem is not None:
                        partit.remove(delElem)
                    new_node = et.Element("endPtTable")
                    new_node.text = re.sub(r'^partitions/', '', finalFileName)
                    place = _get_tag_pos(partit, "endPtTable", ["volume", "endPtTable", "data"], 9999999)
                    partit.insert(place, new_node)

                    header = wellLines[0]

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
                                new_node = et.Element("data")
                                place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                                partit.insert(place, new_node)
                                data = new_node
                                new_node = et.Element("tar", id=targetName)
                                place = _get_tag_pos(data, "tar", partElemLS, 9999999)
                                data.insert(place, new_node)
                            delElem = _get_first_child(partit, "pos")
                            if delElem is not None:
                                data.remove(delElem)
                            new_node = et.Element("pos")
                            new_node.text = str(cPos)
                            place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                            data.insert(place, new_node)
                            delElem = _get_first_child(partit, "neg")
                            if delElem is not None:
                                data.remove(delElem)
                            new_node = et.Element("neg")
                            new_node.text = str(cNeg)
                            place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                            data.insert(place, new_node)
                            delElem = _get_first_child(partit, "undef")
                            if delElem is not None:
                                data.remove(delElem)
                            if cExcl > 0:
                                new_node = et.Element("undef")
                                new_node.text = str(cUndef)
                                place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                                data.insert(place, new_node)
                            delElem = _get_first_child(partit, "excl")
                            if delElem is not None:
                                data.remove(delElem)
                            if cExcl > 0:
                                new_node = et.Element("excl")
                                new_node.text = str(cExcl)
                                place = _get_tag_pos(data, "pos", partElemLS, 9999999)
                                data.insert(place, new_node)

                elif fileformat == "Bio-Rad":
                    wellLines = list(csv.reader(wellfile, delimiter=','))
                    ch1Pos = "0"
                    ch1Neg = "0"
                    ch1sum = 0
                    ch2Pos = "0"
                    ch2Neg = "0"
                    ch2sum = 0

                    if well in headerLookup:
                        if "Ch1" in headerLookup[well] and 1 not in ignoreList:
                            keepCh1 = True
                            header += headerLookup[well]["Ch1"] + "\t" + headerLookup[well]["Ch1"] + "\t"
                        if "Ch2" in headerLookup[well] and 2 not in ignoreList:
                            keepCh2 = True
                            header += headerLookup[well]["Ch2"] + "\t" + headerLookup[well]["Ch2"] + "\t"
                        outTabFile += re.sub(r'\t$', '\n', header)
                    else:
                        headerLookup[well] = {}
                        dyes = ["Ch1", "Ch2"]
                        if len(wellLines) > 1:
                            ch1Pos = ""
                            ch1Neg = ""
                            ch2Pos = ""
                            ch2Neg = ""
                            if re.search(r"\d", wellLines[1][0]) and 1 not in ignoreList:
                                keepCh1 = True
                            if len(wellLines[1]) > 1 and re.search(r"\d", wellLines[1][1]) and 2 not in ignoreList:
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
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh1 = new_node
                            new_node = et.Element("tar", id=headerLookup[well]["Ch1"])
                            place = _get_tag_pos(dataCh1, "tar", partElemLS, 9999999)
                            dataCh1.insert(place, new_node)
                            ch1Pos = ""
                            ch1Neg = ""
                            ch1sum = 2
                        if dataCh2 is None and keepCh2:
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh2 = new_node
                            new_node = et.Element("tar", id=headerLookup[well]["Ch2"])
                            place = _get_tag_pos(dataCh2, "tar", partElemLS, 9999999)
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
                            for splitLine in wellLines[1:]:
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
                                new_node = et.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "pos", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = et.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "neg", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = et.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh1, "neg", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)
                            if keepCh2:
                                new_node = et.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "pos", partElemLS,
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = et.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "neg", partElemLS,
                                                     9999999)
                                dataCh2.insert(place, new_node)

                                new_node = et.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh2, "neg", partElemLS,
                                                     9999999)
                                dataCh2.insert(place, new_node)
                        else:
                            ch1Arr = []
                            ch2Arr = []
                            ch1Cut = 0
                            ch2Cut = 0
                            for splitLine in wellLines[1:]:
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

                            for splitLine in wellLines[1:]:
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
                        new_node = et.Element("endPtTable")
                        new_node.text = re.sub(r'^partitions/', '', finalFileName)
                        place = _get_tag_pos(partit, "endPtTable", ["volume", "endPtTable", "data"], 9999999)
                        partit.insert(place, new_node)
                    else:
                        react.remove(partit)
                elif fileformat == "Stilla":
                    wellLines = list(csv.reader(wellfile, delimiter=','))
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
                        if "Ch1" in headerLookup[well] and 1 not in ignoreList:
                            keepCh1 = True
                            header += headerLookup[well]["Ch1"] + "\t" + headerLookup[well]["Ch1"] + "\t"
                        if "Ch2" in headerLookup[well] and 2 not in ignoreList:
                            keepCh2 = True
                            header += headerLookup[well]["Ch2"] + "\t" + headerLookup[well]["Ch2"] + "\t"
                        if "Ch3" in headerLookup[well] and 3 not in ignoreList:
                            keepCh3 = True
                            header += headerLookup[well]["Ch3"] + "\t" + headerLookup[well]["Ch3"] + "\t"
                        outTabFile += re.sub(r'\t$', '\n', header)
                    else:
                        headerLookup[well] = {}
                        dyes = ["Ch1", "Ch2", "Ch3"]
                        if len(wellLines) > 1:
                            ch1Pos = ""
                            ch1Neg = ""
                            ch2Pos = ""
                            ch2Neg = ""
                            ch3Pos = ""
                            ch3Neg = ""
                            if re.search(r"\d", wellLines[1][0]) and 1 not in ignoreList:
                                keepCh1 = True
                            if len(wellLines[1]) > 1 and re.search(r"\d", wellLines[1][1]) and 2 not in ignoreList:
                                keepCh2 = True
                            if len(wellLines[1]) > 2 and re.search(r"\d", wellLines[1][2]) and 3 not in ignoreList:
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
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh1 = new_node
                            new_node = et.Element("tar", id=headerLookup[well]["Ch1"])
                            place = _get_tag_pos(dataCh1, "tar", partElemLS,
                                                 9999999)
                            dataCh1.insert(place, new_node)
                            ch1Pos = ""
                            ch1Neg = ""
                            ch1sum = 2
                        if dataCh2 is None and keepCh2:
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh2 = new_node
                            new_node = et.Element("tar", id=headerLookup[well]["Ch2"])
                            place = _get_tag_pos(dataCh2, "tar", partElemLS,
                                                 9999999)
                            dataCh2.insert(place, new_node)
                            ch2Pos = ""
                            ch2Neg = ""
                            ch2sum = 2
                        if dataCh3 is None and keepCh3:
                            new_node = et.Element("data")
                            place = _get_tag_pos(partit, "data", ["volume", "endPtTable", "data"], 9999999)
                            partit.insert(place, new_node)
                            dataCh3 = new_node
                            new_node = et.Element("tar", id=headerLookup[well]["Ch3"])
                            place = _get_tag_pos(dataCh3, "tar", partElemLS,
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
                            for splitLine in wellLines[1:]:
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
                                new_node = et.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "pos", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = et.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh1, "neg", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)

                                new_node = et.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh1, "neg", partElemLS, 9999999)
                                dataCh1.insert(place, new_node)
                            if keepCh2:
                                new_node = et.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "pos", partElemLS, 9999999)
                                dataCh2.insert(place, new_node)

                                new_node = et.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh2, "neg", partElemLS, 9999999)
                                dataCh2.insert(place, new_node)

                                new_node = et.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh2, "neg", partElemLS, 9999999)
                                dataCh2.insert(place, new_node)
                            if keepCh3:
                                new_node = et.Element("pos")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh3, "pos", partElemLS, 9999999)
                                dataCh3.insert(place, new_node)

                                new_node = et.Element("neg")
                                new_node.text = "0"
                                place = _get_tag_pos(dataCh3, "neg", partElemLS, 9999999)
                                dataCh3.insert(place, new_node)

                                new_node = et.Element("undef")
                                new_node.text = str(countPart)
                                place = _get_tag_pos(dataCh3, "neg", partElemLS, 9999999)
                                dataCh3.insert(place, new_node)
                        else:
                            ch1Arr = []
                            ch2Arr = []
                            ch3Arr = []
                            ch1Cut = 0
                            ch2Cut = 0
                            ch3Cut = 0
                            for splitLine in wellLines[1:]:
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

                            for splitLine in wellLines[1:]:
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
                        new_node = et.Element("endPtTable")
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
        tarTypeLookup = {}
        tarDyeLookup = {}

        transSamTar = _sampleTypeToDics(rootEl._node)
        targets = _get_all_children(rootEl._node, "target")
        for target in targets:
            if target.attrib['id'] != "":
                tarId = target.attrib['id']
                forType = _get_first_child_text(target, "type")
                if forType != "":
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
            pFileName = ""
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    pSample = forId.attrib['id']
            partit = _get_first_child(react, "partitions")
            if partit is not None:
                endPtTable = _get_first_child_text(partit, "endPtTable")
                if endPtTable != "":
                    pFileName = endPtTable
                pVolume = _get_first_child_text(partit, "volume")
                partit_datas = _get_all_children(partit, "data")
                for partit_data in partit_datas:
                    pSampleType = ""
                    pTarget = ""
                    pTargetType = ""
                    pDye = ""
                    forId = _get_first_child(partit_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] != "":
                            pTarget = forId.attrib['id']
                            pSampleType = transSamTar[pSample][pTarget]
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
        retVal = ""

        # Get the position number if required
        wellPos = str(reactPos)
        if re.search(r"\D\d+", wellPos):
            old_letter = ord(re.sub(r"\d", "", wellPos.upper())) - ord("A")
            old_nr = int(re.sub(r"\D", "", wellPos))
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

    def setReactionVolume(self, vol):
        """Imports data from a tab seperated table file with react fluorescence data.

        Args:
            self: The class self parameter.
            vol: The default volume of each reaction in microliter.

        Returns:
            None.
        """

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        rdml_version = rootPar.get('version')
        if rdml_version not in ['1.4']:
            return

        if vol == "":
            return

        vol = re.sub(r"[^0-9\.]", "", vol)
        try:
            mVol = float(vol)
        except ValueError:
            return
        if not math.isfinite(mVol):
            return
        if mVol <= 0.0:
            return

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            _change_subelement(react, "vol", ["sample", "vol", "data", "partitions"], vol, True, "float")

    def getreactjson(self, curves=True):
        """Returns a json of the react data including fluorescence data.

        Args:
            self: The class self parameter.
            curves: Include amplification and melting curves

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
        anyCorrections = 0

        for react in reacts:
            react_json = {
                "id": react.get('id'),
            }
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    react_json["sample"] = forId.attrib['id']
            react_json["vol"] = _get_first_child_text(react, "vol")
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
                _add_first_child_to_dic(react_data, in_react, True, "N0")
                _add_first_child_to_dic(react_data, in_react, True, "Ncopy")
                _add_first_child_to_dic(react_data, in_react, True, "ampEffMet")
                _add_first_child_to_dic(react_data, in_react, True, "ampEff")
                _add_first_child_to_dic(react_data, in_react, True, "ampEffSE")
                _add_first_child_to_dic(react_data, in_react, True, "corrF")
                _add_first_child_to_dic(react_data, in_react, True, "corrP")
                _add_first_child_to_dic(react_data, in_react, True, "bgFluor")
                # Calculate the correction factors
                corrFac = _get_first_child_text(react_data, "corrF")
                calcCorr = 1.0
                if not corrFac == "":
                    try:
                        calcCorr = float(corrFac)
                    except ValueError:
                        calcCorr = 1.0
                    if not math.isfinite(calcCorr):
                        calcCorr = 1.0
                    if calcCorr > 1.0:
                        calcCorr = 1.0
                plateFac = _get_first_child_text(react_data, "corrP")
                calcPlate = 1.0
                if not plateFac == "":
                    try:
                        calcPlate = float(plateFac)
                    except ValueError:
                        calcPlate = 0.0
                    if not math.isfinite(calcPlate):
                        calcCorr = 0.0
                    if calcPlate == 0.0:
                        calcCorr = 0.0
                    calcCorr = calcCorr / calcPlate
                calcN0 = _get_first_child_text(react_data, "N0")
                in_react["corrN0"] = calcN0
                if not calcN0 == "":
                    try:
                        calcN0 = float(calcN0)
                    except ValueError:
                        pass
                    else:
                        if math.isfinite(calcN0):
                            if calcCorr > 0.0001:
                                finalN0 = calcCorr * calcN0
                                in_react["corrN0"] = "{:.2e}".format(finalN0)
                                if not corrFac == "":
                                    anyCorrections = 1
                                if not plateFac == "":
                                    anyCorrections = 1
                calcNcopy = _get_first_child_text(react_data, "Ncopy")
                in_react["corrNcopy"] = calcNcopy
                if not calcNcopy == "":
                    try:
                        calcNcopy = float(calcNcopy)
                    except ValueError:
                        pass
                    else:
                        if math.isfinite(calcNcopy):
                            if calcCorr > 0.0001:
                                finalNcopy = calcCorr * calcNcopy
                                in_react["corrNcopy"] = "{:.4f}".format(finalNcopy)
                                if not corrFac == "":
                                    anyCorrections = 1
                                if not plateFac == "":
                                    anyCorrections = 1
                if calcCorr <= 0.0001:
                    if calcCorr < 0.0:
                        in_react["corrN0"] = "-1.0"
                        in_react["corrNcopy"] = "-1.0"
                    else:
                        in_react["corrN0"] = "0.0"
                        in_react["corrNcopy"] = "0.0"
                if calcPlate < 0.0:
                    in_react["corrN0"] = "-1.0"
                    in_react["corrNcopy"] = "-1.0"

                _add_first_child_to_dic(react_data, in_react, True, "meltTemp")
                _add_first_child_to_dic(react_data, in_react, True, "excl")
                _add_first_child_to_dic(react_data, in_react, True, "note")
                _add_first_child_to_dic(react_data, in_react, True, "endPt")
                _add_first_child_to_dic(react_data, in_react, True, "bgFluor")
                _add_first_child_to_dic(react_data, in_react, True, "bgFluorSlp")
                _add_first_child_to_dic(react_data, in_react, True, "quantFluor")
                if curves:
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
                if endPtTable != "":
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
                    _add_first_child_to_dic(partit_data, in_partit, True, "excluded")
                    _add_first_child_to_dic(partit_data, in_partit, True, "note")
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
        all_data["anyCalcCorrections"] = anyCorrections
        all_data["adp_fluor_min"] = adp_fluor_min
        all_data["adp_fluor_max"] = adp_fluor_max
        all_data["mdp_tmp_min"] = mdp_tmp_min
        all_data["mdp_tmp_max"] = mdp_tmp_max
        all_data["mdp_fluor_min"] = mdp_fluor_min
        all_data["mdp_fluor_max"] = mdp_fluor_max
        all_data["max_data_len"] = max_data
        all_data["max_partition_data_len"] = max_partition_data
        return all_data

    def setVolume(self, vReact, vol):
        """Adds a volume to the reaction from the RDML data.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vol: The volume in microliters.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                _change_subelement(react, "vol", ["sample", "vol", "data", "partitions"], vol, True, "float")
                break
        return

    def setVolumeGrp(self, vWells, vol):
        """Adds a volume to several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vol: The volume in microliters.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setVolume(well, vol)
        return

    def removeReact(self, vReact):
        """Removes the reaction from the RDML data.

        Args:
            self: The class self parameter.
            vReact: The reaction id.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                self._node.remove(react)
                break
        return

    def removeReactGrp(self, vWells):
        """Removes several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.removeReact(well)
        return

    def removeClasReactTar(self, vReact, vTar):
        """Removes the tar data of a reaction from the RDML data.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                remain_data = 0
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    forId = _get_first_child(react_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] == vTar:
                            react.remove(react_data)
                            break
                remain_data += len(_get_all_children(react, "data"))
                partit = _get_first_child(react, "partitions")
                if partit is not None:
                    remain_data += len(_get_all_children(partit, "data"))
                if remain_data == 0:
                    self._node.remove(react)
                    break
        return

    def removeClasReactTarGrp(self, vWells, vTar):
        """Removes the tar data of several reactions from the RDML data.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.removeClasReactTar(well, vTar)
        return

    def removeDigiReactTar(self, vReact, vTar):
        """Removes the tar data of a reaction from the RDML data.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                remain_data = 0
                partit = _get_first_child(react, "partitions")
                if partit is not None:
                    part_datas = _get_all_children(partit, "data")
                    for part_data in part_datas:
                        forId = _get_first_child(part_data, "tar")
                        if forId is not None:
                            if forId.attrib['id'] == vTar:
                                partit.remove(part_data)
                                break
                remain_data += len(_get_all_children(react, "data"))
                partit = _get_first_child(react, "partitions")
                if partit is not None:
                    remain_data += len(_get_all_children(partit, "data"))
                if remain_data == 0:
                    self._node.remove(react)
                    break
        return

    def removeDigiReactTarGrp(self, vWells, vTar):
        """Removes the tar data of several reactions from the RDML data.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.removeDigiReactTar(well, vTar)
        return

    def setClasSample(self, vReact, vNewSam):
        """Changes the sample id for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vNewTar: The new target id to change to.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                forId = _get_first_child(react, "sample")
                if forId is not None:
                    forId.attrib['id'] = vNewSam
        return

    def setClasSampleGrp(self, vWells, vNewSam):
        """Changes the target id for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vNewTar: The new target id to change to.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setClasSample(well, vNewSam)
        return

    def setClasTarget(self, vReact, vTar, vNewTar):
        """Changes the target id for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.
            vNewTar: The new target id to change to.

        Returns:
            Nothing, updates RDML data.
        """

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    forId = _get_first_child(react_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] == vTar:
                            forId.attrib['id'] = vNewTar
        return

    def setClasTragetGrp(self, vWells, vTar, vNewTar):
        """Changes the target id for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.
            vNewTar: The new target id to change to.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setClasTarget(well, vTar, vNewTar)
        return

    def setClasExcl(self, vReact, vTar, vExcl, vAppend):
        """Saves the excl string for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.
            vExcl: The exclusion string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        dataXMLelements = _getXMLDataType()
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    forId = _get_first_child(react_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] == vTar:
                            oldData = ""
                            if vAppend:
                                oldData = _get_first_child_text(react_data, "excl")
                            _change_subelement(react_data, "excl", dataXMLelements, oldData + vExcl, True, "string")
        return

    def setClasExclGrp(self, vWells, vTar, vExcl, vAppend):
        """Saves the excl string for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.
            vExcl: The exclusion string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setClasExcl(well, vTar, vExcl, vAppend)
        return

    def setDigiExcl(self, vReact, vTar, vExcl, vAppend):
        """Saves the excl string for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.
            vExcl: The exclusion string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        ver = rootPar.get('version')
        if ver in ["1.0", "1.1", "1.2"]:
            return

        dataXMLelements = _getXMLPartitionDataType()
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                partit = _get_first_child(react, "partitions")
                if partit is not None:
                    part_datas = _get_all_children(partit, "data")
                    for part_data in part_datas:
                        forId = _get_first_child(part_data, "tar")
                        if forId is not None:
                            if forId.attrib['id'] == vTar:
                                oldData = ""
                                if vAppend:
                                    oldData = _get_first_child_text(part_data, "excluded")
                                _change_subelement(part_data, "excluded", dataXMLelements, oldData + vExcl, True, "string")
        return

    def setDigiExclGrp(self, vWells, vTar, vExcl, vAppend):
        """Saves the excl string for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.
            vExcl: The exclusion string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setDigiExcl(well, vTar, vExcl, vAppend)
        return

    def setClasNote(self, vReact, vTar, vNote, vAppend):
        """Saves the note string for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.
            vNote: The note string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        ver = rootPar.get('version')
        dataXMLelements = _getXMLDataType()

        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                react_datas = _get_all_children(react, "data")
                for react_data in react_datas:
                    forId = _get_first_child(react_data, "tar")
                    if forId is not None:
                        if forId.attrib['id'] == vTar:
                            if ver in ["1.3", "1.4"]:
                                oldData = ""
                                if vAppend:
                                    oldData = _get_first_child_text(react_data, "note")
                                _change_subelement(react_data, "note", dataXMLelements, oldData + vNote, True, "string")
        return

    def setClasNoteGrp(self, vWells, vTar, vNote, vAppend):
        """Saves the note string for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.
            vNote: The note string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setClasNote(well, vTar, vNote, vAppend)
        return

    def setDigiNote(self, vReact, vTar, vNote, vAppend):
        """Saves the note string for one react/data combination.

        Args:
            self: The class self parameter.
            vReact: The reaction id.
            vTar: The target id.
            vNote: The note string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        ver = rootPar.get('version')
        if ver in ["1.0", "1.1", "1.2"]:
            return

        dataXMLelements = _getXMLPartitionDataType()
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            if int(react.get('id')) == int(vReact):
                partit = _get_first_child(react, "partitions")
                if partit is not None:
                    part_datas = _get_all_children(partit, "data")
                    for part_data in part_datas:
                        forId = _get_first_child(part_data, "tar")
                        if forId is not None:
                            if forId.attrib['id'] == vTar:
                                oldData = ""
                                if vAppend:
                                    oldData = _get_first_child_text(part_data, "note")
                                _change_subelement(part_data, "note", dataXMLelements, oldData + vNote, True, "string")
        return

    def setDigiNoteGrp(self, vWells, vTar, vNote, vAppend):
        """Saves the note string for several react/data combination.

        Args:
            self: The class self parameter.
            vWells: The a range of wells like "B2-C4".
            vTar: The target id.
            vNote: The note string to save.
            vAppend: If True, append to excisting string, if False replace.

        Returns:
            Nothing, updates RDML data.
        """

        wells = _wellRangeToList(vWells, self["pcrFormat_columns"])
        for well in wells:
            self.setDigiNote(well, vTar, vNote, vAppend)
        return

    def webAppLinRegPCR(self, pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=False, excludeNoPlateau=True,
                        excludeEfficiency="outlier", excludeInstableBaseline=True):
        """Performs LinRegPCR on the run. Modifies the cq values and returns a json with additional data.

        Args:
            self: The class self parameter.
            pcrEfficiencyExl: Exclude samples with an efficiency outside the given range (0.05).
            updateTargetEfficiency: If true, the PCR efficiencies and errors of the targets are updated.
            updateRDML: If true, update the RDML data with the calculated values.
            excludeNoPlateau: If true, samples without plateau are excluded from mean PCR efficiency calculation.
            excludeEfficiency: Choose "outlier", "mean", "include" to exclude based on indiv PCR eff.
            excludeInstableBaseline: If true, instable baselines give a baseline error.

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            rawData: A 2d array with the raw fluorescence values
            baselineCorrectedData: A 2d array with the baseline corrected raw fluorescence values
            resultsList: A 2d array object.
            resultsCSV: A csv string.
        """

        allData = self.getreactjson()
        res = self.linRegPCR(pcrEfficiencyExl=pcrEfficiencyExl,
                             updateTargetEfficiency=updateTargetEfficiency,
                             updateRDML=updateRDML,
                             excludeNoPlateau=excludeNoPlateau,
                             excludeEfficiency=excludeEfficiency,
                             excludeInstableBaseline=excludeInstableBaseline,
                             saveRaw=False,
                             saveBaslineCorr=True,
                             saveResultsList=True,
                             saveResultsCSV=False,
                             verbose=False)
        if "baselineCorrectedData" in res:
            bas_cyc_max = len(res["baselineCorrectedData"][0]) - 5
            bas_fluor_min = 99999999
            bas_fluor_max = 0.0
            for row in range(1, len(res["baselineCorrectedData"])):
                bass_json = []
                for col in range(5, len(res["baselineCorrectedData"][row])):
                    cyc = res["baselineCorrectedData"][0][col]
                    fluor = res["baselineCorrectedData"][row][col]
                    if np.isfinite(fluor) and fluor > 0.0:
                        bas_fluor_min = min(bas_fluor_min, float(fluor))
                        bas_fluor_max = max(bas_fluor_max, float(fluor))
                        in_bas = [cyc, fluor, ""]
                        bass_json.append(in_bas)
                # Fixme do not loop over all, use sorted data and clever moving
                for react in allData["reacts"]:
                    if react["id"] == res["baselineCorrectedData"][row][0]:
                        for data in react["datas"]:
                            if data["tar"] == res["baselineCorrectedData"][row][3]:
                                data["bass"] = list(bass_json)
            allData["bas_cyc_max"] = bas_cyc_max
            allData["bas_fluor_min"] = bas_fluor_min
            allData["bas_fluor_max"] = bas_fluor_max

        if "resultsList" in res:
            header = res["resultsList"].pop(0)
            resList = sorted(res["resultsList"], key=_sort_list_int)
            for rRow in range(0, len(resList)):
                for rCol in range(0, len(resList[rRow])):
                    if isinstance(resList[rRow][rCol], np.float64) and not np.isfinite(resList[rRow][rCol]):
                        resList[rRow][rCol] = ""
                    if isinstance(resList[rRow][rCol], float) and not math.isfinite(resList[rRow][rCol]):
                        resList[rRow][rCol] = ""
            allData["LinRegPCR_Result_Table"] = json.dumps([header] + resList, cls=NpEncoder)

        if "error" in res:
            allData["error"] = res["error"]

        if "resultsLinRegPCRReport" in res:
            allData["resultsLinRegPCRReport"] = res["resultsLinRegPCRReport"]

        return allData

    def linRegPCR(self, pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=False,
                  excludeNoPlateau=True, excludeEfficiency="outlier",
                  excludeInstableBaseline=True, commaConv=False, ignoreExclusion=False,
                  saveRaw=False, saveBaslineCorr=False, saveResultsList=False, saveResultsCSV=False,
                  timeRun=False, verbose=False):
        """Performs LinRegPCR on the run. Modifies the cq values and returns a json with additional data.

        Args:
            self: The class self parameter.
            pcrEfficiencyExl: Exclude samples with an efficiency outside the given range (0.05).
            updateTargetEfficiency: If true, the PCR efficiencies and errors of the targets are updated.
            updateRDML: If true, update the RDML data with the calculated values.
            excludeNoPlateau: If true, samples without plateau are excluded from mean PCR efficiency calculation.
            excludeEfficiency: Choose "outlier", "mean", "include" to exclude based on indiv PCR eff.
            excludeInstableBaseline: If true, instable baselines give a baseline error.
            commaConv: If true, convert comma separator to dot.
            ignoreExclusion: If true, ignore the RDML exclusion strings.
            saveRaw: If true, no raw values are given in the returned data
            saveBaslineCorr: If true, no baseline corrected values are given in the returned data
            saveResultsList: If true, return a 2d array object.
            saveResultsCSV: If true, return a csv string.
            timeRun: If true, print runtime for baseline and total.
            verbose: If true, comment every performed step.

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            rawData: A 2d array with the raw fluorescence values
            baselineCorrectedData: A 2d array with the baseline corrected raw fluorescence values
            resultsList: A 2d array object.
            resultsCSV: A csv string.
        """

        try:
            pcrEfficiencyExl = float(pcrEfficiencyExl)
        except ValueError:
            raise RdmlError('LinRegPCR pcrEfficiencyExl parameter must be a float number.')
        if not np.isfinite(pcrEfficiencyExl):
            raise RdmlError('LinRegPCR pcrEfficiencyExl parameter must be a float number.')
        if excludeEfficiency not in ["outlier", "mean", "include"]:
            raise RdmlError('LinRegPCR excludeEfficiency parameter must be "outlier", "mean" or "include".')

        if type(updateTargetEfficiency) != bool:
            raise RdmlError('LinRegPCR updateTargetEfficiency parameter must be True or False.')
        if type(updateRDML) != bool:
            raise RdmlError('LinRegPCR updateRDML parameter must be True or False.')
        if type(excludeNoPlateau) != bool:
            raise RdmlError('LinRegPCR excludeNoPlateau parameter must be True or False.')
        if type(excludeInstableBaseline) != bool:
            raise RdmlError('LinRegPCR excludeInstableBaseline parameter must be True or False.')
        if type(commaConv) != bool:
            raise RdmlError('LinRegPCR commaConv parameter must be True or False.')
        if type(ignoreExclusion) != bool:
            raise RdmlError('LinRegPCR ignoreExclusion parameter must be True or False.')
        if type(saveRaw) != bool:
            raise RdmlError('LinRegPCR saveRaw parameter must be True or False.')
        if type(saveBaslineCorr) != bool:
            raise RdmlError('LinRegPCR saveBaslineCorr parameter must be True or False.')
        if type(saveResultsList) != bool:
            raise RdmlError('LinRegPCR saveResultsList parameter must be True or False.')
        if type(saveResultsCSV) != bool:
            raise RdmlError('LinRegPCR saveResultsCSV parameter must be True or False.')
        if type(timeRun) != bool:
            raise RdmlError('LinRegPCR timeRun parameter must be True or False.')
        if type(verbose) != bool:
            raise RdmlError('LinRegPCR verbose parameter must be True or False.')

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        dataVersion = rootPar.get('version')
        transSamTar = _sampleTypeToDics(rootPar)

        if dataVersion not in ["1.4"]:
            raise RdmlError('LinRegPCR requires RDML version >= 1.4. Got v' + dataVersion)
        err = ""

        ##############################
        # Collect the data in arrays #
        ##############################

        # res is a 2 dimensional array accessed only by
        # variables, so columns might be added here
        header = [["id",  # 0
                   "well",  # 1
                   "sample",  # 2
                   "sample type",  # 3
                   "sample nucleotide",  # 4
                   "target",   # 5
                   "target chemistry",  # 6
                   "excluded",   # 7
                   "note",   # 8
                   "baseline",   # 9
                   "plateau",   # 10
                   "plateau / baseline",   # 11
                   "lower limit",   # 12
                   "upper limit",   # 13
                   "n in log phase",   # 14
                   "last log cycle",   # 15
                   "n included",   # 16
                   "indiv PCR eff",   # 17
                   "R2",   # 18
                   "PCR eff",   # 19
                   "standard error of PCR eff",   # 20
                   "threshold",  # 21
                   "Cq",   # 22
                   "indiv Ncopy",   # 23
                   "Ncopy",   # 24
                   "amplification",   # 25
                   "baseline error",   # 26
                   "instable baseline",   # 27
                   "plateau",   # 28
                   "noisy sample",   # 29
                   "excl PCR efficiency",   # 30
                   "short log lin phase",   # 31
                   "Cq is shifting",   # 32
                   "too low Cq eff",   # 33
                   "too low Cq N0",   # 34
                   "used for W-o-L setting",   # 35
                   "nAmpli",   # 36
				   "vol",   # 37
				   "dNTPs",   # 38
				   "dyeConc",   # 39
				   "primer len for",   # 40
				   "primer conc for",   # 41
				   "primer len rev",   # 42
				   "primer conc rev",   # 43
				   "probe1 len",   # 44
				   "probe1 conc",   # 45
				   "probe2 len",   # 46
				   "probe2 conc",   # 47
				   "amplicon len",
                   "del threshold",
                   "del N0",
                   "del ncopy",
                   "del fact"
                   ]]   # 48
        rar_id = 0
        rar_well = 1
        rar_sample = 2
        rar_sample_type = 3
        rar_sample_nucleotide = 4
        rar_tar = 5
        rar_tar_chemistry = 6
        rar_excl = 7
        rar_note = 8
        rar_baseline = 9
        rar_plat_val = 10
        rar_plat_base = 11
        rar_lower_limit = 12
        rar_upper_limit = 13
        rar_n_log = 14
        rar_stop_log = 15
        rar_n_included = 16
        rar_indiv_PCR_eff = 17
        rar_R2 = 18
        rar_PCR_eff = 19
        rar_PCR_eff_err = 20
        rar_threshold = 21
        rar_Cq = 22
        rar_indiv_Ncopy = 23
        rar_Ncopy = 24
        rar_amplification = 25
        rar_baseline_error = 26
        rar_instable_baseline = 27
        rar_plateau = 28
        rar_noisy_sample = 29
        rar_excl_eff = 30
        rar_shortLogLinPhase = 31
        rar_CqIsShifting = 32
        rar_tooLowCqEff = 33
        rar_tooLowCqN0 = 34
        rar_isUsedInWoL = 35
        rar_nAmpli = 36
        rar_vol = 37
        rar_dNTPs = 38
        rar_dyeConc = 39
        rar_primer_len_for = 40
        rar_primer_conc_for = 41
        rar_primer_len_rev = 42
        rar_primer_conc_rev = 43
        rar_probe1_len = 44
        rar_probe1_conc = 45
        rar_probe2_len = 46
        rar_probe2_conc = 47
        rar_amplicon_len = 48
        del_threshold = 49
        del_N0 = 50
        del_ncopy = 51
        del_fact = 52

        copyThreshold = _get_copies_threshold()

        # Constand default values
        AVOGADRO = 6.02214076e23
        DEF_VOL_96 = 20.0
        DEF_VOL_384 = 10.0
        DEF_DNTP = 200.0
        DEF_DYE_CONC = 0.98
        DEF_PRIMER_CONC = 250.0
        DEF_PRIMER_LEN = 20.0
        DEF_AMPLICON_LEN = 100.0

        #############################################
        #  Read all data from RDML into the arrays  #
        #############################################

        res = []
        finalData = {}

        # First get the max number of cycles and create the numpy array
        colCount = 0
        adp_cyc_min = 1000000000.0
        adp_cyc_max = 0.0
        reacts = _get_all_children(self._node, "react")
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                colCount += 1
                adps = _get_all_children(react_data, "adp")
                for adp in adps:
                    cyc = float(_get_first_child_text(adp, "cyc"))
                    adp_cyc_min = min(adp_cyc_min, cyc)
                    adp_cyc_max = max(adp_cyc_max, cyc)
        # This spanns the matrix from min to max (len = min + max + 1)
        # Gaps in this matrix will be filled or cause an error
        # The matrix is guranteed to be with valid values
        adp_cyc_min = int(math.ceil(adp_cyc_min))
        adp_cyc_max = int(math.ceil(adp_cyc_max))
        if adp_cyc_min < 1:
            raise RdmlError('adp cycles must be >= 1.')

        # print(str(adp_cyc_min) + " - " + str(adp_cyc_max))

        # spFl is the shape for all fluorescence numpy data arrays
        # cycle number is arrar pos + 1 - zero is first
        spFl = (colCount, adp_cyc_max - adp_cyc_min + 1)
        spNum = np.arange(0, colCount)
        # print(spFl)
        rawFluor = np.zeros(spFl, dtype=np.float64)
        rawFluor[rawFluor <= 0.00000001] = np.nan

        # Initialization of the vecNoAmplification vector
        vecExcludedByUser = np.zeros(spFl[0], dtype=np.bool_)
        rdmlElemData = []

        # Now process the data for numpy and create results array
        anyRawData = False
        rowCount = 0
        reactVol = np.zeros(spFl[0], dtype=np.float64)
        volSum = 0.0
        volNum = 0
        devVol = DEF_VOL_384
        if int(self["pcrFormat_columns"]) < 13:
            devVol = DEF_VOL_96
        for react in reacts:
            posId = react.get('id')
            pIdNumber = (int(posId) - 1) % int(self["pcrFormat_columns"]) + 1
            pIdLetter = chr(ord("A") + int((int(posId) - 1) / int(self["pcrFormat_columns"])))
            pWell = pIdLetter + str(pIdNumber)
            sample = ""
            forId = _get_first_child(react, "sample")
            if forId is not None:
                if forId.attrib['id'] != "":
                    sample = forId.attrib['id']
            volTemp = _save_float(_get_first_child_text(react, "vol"))
            if volTemp != "":
                volSum += volTemp
                volNum += 1
                reactVol[rowCount] = volTemp
            else:
                reactVol[rowCount] = devVol
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                forId = _get_first_child(react_data, "tar")
                target = ""
                if forId is not None:
                    if forId.attrib['id'] != "":
                        target = forId.attrib['id']

                if ignoreExclusion:
                    excl = ""
                else:
                    excl = _get_first_child_text(react_data, "excl")
                    excl = _cleanErrorString(excl, "amp")
                    excl = re.sub(r'^;|;$', '', excl)
                if not excl == "":
                    vecExcludedByUser[rowCount] = True

                noteVal = _get_first_child_text(react_data, "note")
                noteVal = _cleanErrorString(noteVal, "amp")
                noteVal = re.sub(r'^;|;$', '', noteVal)
                # To update the RDML data later
                rdmlElemData.append(react_data)
                res.append([posId, pWell, sample, "",  "",  target, "", excl, noteVal, "",
                            "", "", "", "", "",  "", "", "", "", "",
                            "", "", "", "", "",  "", "", "", "", "",
                            "", "", "", "", "",  "", "", "", "", "",
                            "", "", "", "", "",  "", "", "", "", "", "", "", "" ])  # Must match header length
                adps = _get_all_children(react_data, "adp")
                for adp in adps:
                    cyc = int(math.ceil(float(_get_first_child_text(adp, "cyc"))))
                    fluor = _get_first_child_text(adp, "fluor")
                    if commaConv:
                        noDot = fluor.replace(".", "")
                        fluor = noDot.replace(",", ".")
                    rawFluor[rowCount, cyc - adp_cyc_min] = float(fluor)
                    anyRawData = True
                rowCount += 1
        averageVolume = devVol
        if volNum > 0:
            averageVolume = volSum / float(volNum)
        if anyRawData == False:
            raise RdmlError('LinRegPCR requires raw data. No raw data were found in this run.')

        # Look up sample and target information
        parExp = self._node.getparent()
        parRoot = parExp.getparent()

        dicLU_dyes = {}
        dicLU_dye_conc = {}
        dicLU_dye_dNTPs = {}
        luDyes = _get_all_children(parRoot, "dye")
        for lu_dye in luDyes:
            lu_chemistry = _get_first_child_text(lu_dye, "dyeChemistry")
            if lu_chemistry == "":
                lu_chemistry = "non-saturating DNA binding dye"
            if lu_dye.attrib['id'] != "":
                dicLU_dyes[lu_dye.attrib['id']] = lu_chemistry
            dicLU_dye_conc[lu_dye.attrib['id']] = DEF_DYE_CONC
            lu_dyeConc = _get_first_child_text(lu_dye, "dyeConc")
            dyeConcVal = _save_float(lu_dyeConc)
            if dyeConcVal != "":
                dicLU_dye_conc[lu_dye.attrib['id']] = dyeConcVal
            dicLU_dye_dNTPs[lu_dye.attrib['id']] = DEF_DNTP
            lu_dNTPs = _get_first_child_text(lu_dye, "dyeConc")
            dyeConcDNTPs = _save_float(lu_dNTPs)
            if dyeConcDNTPs != "":
                dicLU_dye_dNTPs[lu_dye.attrib['id']] = dyeConcDNTPs
        dicLU_targets = {}
        dicLU_target_dyeConc = {}
        dicLU_target_dNTPs = {}
        dicLU_target_fw_len = {}
        dicLU_target_rv_len = {}
        dicLU_target_probe1_len = {}
        dicLU_target_probe2_len = {}
        dicLU_target_amp_len = {}
        dicLU_target_fw_conc = {}
        dicLU_target_rv_conc = {}
        dicLU_target_probe1_conc = {}
        dicLU_target_probe2_conc = {}
        luTargets = _get_all_children(parRoot, "target")
        for lu_target in luTargets:
            if lu_target.attrib['id'] == "":
                continue
            forId = _get_first_child(lu_target, "dyeId")
            lu_dyeId = ""
            if forId is not None:
                if forId.attrib['id'] != "":
                    lu_dyeId = forId.attrib['id']
            dicLU_targets[lu_target.attrib['id']] = "non-saturating DNA binding dye"
            dicLU_target_dyeConc[lu_target.attrib['id']] = DEF_DYE_CONC
            dicLU_target_dNTPs[lu_target.attrib['id']] = DEF_DNTP
            dicLU_target_fw_len[lu_target.attrib['id']] = DEF_PRIMER_LEN
            dicLU_target_rv_len[lu_target.attrib['id']] = DEF_PRIMER_LEN
            dicLU_target_probe1_len[lu_target.attrib['id']] = DEF_PRIMER_LEN
            dicLU_target_probe2_len[lu_target.attrib['id']] = DEF_PRIMER_LEN
            dicLU_target_amp_len[lu_target.attrib['id']] = DEF_AMPLICON_LEN
            dicLU_target_fw_conc[lu_target.attrib['id']] = DEF_PRIMER_CONC
            dicLU_target_rv_conc[lu_target.attrib['id']] = DEF_PRIMER_CONC
            dicLU_target_probe1_conc[lu_target.attrib['id']] = DEF_PRIMER_CONC
            dicLU_target_probe2_conc[lu_target.attrib['id']] = DEF_PRIMER_CONC
            if lu_dyeId != "" and lu_dyeId in dicLU_dyes:
                dicLU_targets[lu_target.attrib['id']] = dicLU_dyes[lu_dyeId]
                dicLU_target_dyeConc[lu_target.attrib['id']] = dicLU_dye_conc[lu_dyeId]
                dicLU_target_dNTPs[lu_target.attrib['id']] = dicLU_dye_dNTPs[lu_dyeId]
            seqNode = _get_first_child(lu_target, "sequences")
            if seqNode != None:
                fwPrimNode = _get_first_child(seqNode, "forwardPrimer")
                if fwPrimNode != None:
                    fwPrimLenTmp = len(str(_get_first_child_text(fwPrimNode, "sequence")))
                    if fwPrimLenTmp > 5:
                        dicLU_target_fw_len[lu_target.attrib['id']] = float(fwPrimLenTmp)
                    fwPrimConcTmp = _save_float(_get_first_child_text(fwPrimNode, "oligoConc"))
                    if fwPrimConcTmp != "":
                        dicLU_target_fw_conc[lu_target.attrib['id']] = fwPrimConcTmp
                rvPrimNode = _get_first_child(seqNode, "reversePrimer")
                if rvPrimNode != None:
                    rvPrimLenTmp = len(str(_get_first_child_text(rvPrimNode, "sequence")))
                    if rvPrimLenTmp > 5:
                        dicLU_target_rv_len[lu_target.attrib['id']] = float(rvPrimLenTmp)
                    rvPrimConcTmp = _save_float(_get_first_child_text(rvPrimNode, "oligoConc"))
                    if rvPrimConcTmp != "":
                        dicLU_target_rv_conc[lu_target.attrib['id']] = rvPrimConcTmp
                probe1PrimNode = _get_first_child(seqNode, "probe1")
                if probe1PrimNode != None:
                    probe1PrimLenTmp = len(str(_get_first_child_text(probe1PrimNode, "sequence")))
                    if probe1PrimLenTmp > 5:
                        dicLU_target_probe1_len[lu_target.attrib['id']] = float(probe1PrimLenTmp)
                    probe1PrimConcTmp = _save_float(_get_first_child_text(probe1PrimNode, "oligoConc"))
                    if probe1PrimConcTmp != "":
                        dicLU_target_probe1_conc[lu_target.attrib['id']] = probe1PrimConcTmp
                probe2PrimNode = _get_first_child(seqNode, "probe2")
                if probe2PrimNode != None:
                    probe2PrimLenTmp = len(str(_get_first_child_text(probe2PrimNode, "sequence")))
                    if probe2PrimLenTmp > 5:
                        dicLU_target_probe2_len[lu_target.attrib['id']] = float(probe2PrimLenTmp)
                    probe2PrimConcTmp = _save_float(_get_first_child_text(probe2PrimNode, "oligoConc"))
                    if probe2PrimConcTmp != "":
                        dicLU_target_probe2_conc[lu_target.attrib['id']] = probe2PrimConcTmp
                ampPrimNode = _get_first_child(seqNode, "amplicon")
                if ampPrimNode != None:
                    ampPrimLenTmp = len(str(_get_first_child_text(ampPrimNode, "sequence")))
                    if ampPrimLenTmp > 30:
                        dicLU_target_amp_len[lu_target.attrib['id']] = float(ampPrimLenTmp)

        dicLU_samNucl = {}
        luSamples = _get_all_children(parRoot, "sample")
        for lu_sample in luSamples:
            lu_Nucl = _string_to_bool(_get_first_child_text(lu_sample, "doubleStranded"), False)
            if lu_sample.attrib['id'] != "":
                dicLU_samNucl[lu_sample.attrib['id']] = lu_Nucl

        # Update the table with dictionary help
        for oRow in range(0, spFl[0]):
            if res[oRow][rar_sample] != "":
                if res[oRow][rar_tar] != "":
                    res[oRow][rar_sample_type] = transSamTar[res[oRow][rar_sample]][res[oRow][rar_tar]]
                if dicLU_samNucl[res[oRow][rar_sample]] == True:
                    res[oRow][rar_sample_nucleotide] = "ds"
                else:
                    res[oRow][rar_sample_nucleotide] = "ss"
            if res[oRow][rar_tar] != "":
                res[oRow][rar_tar_chemistry] = dicLU_targets[res[oRow][rar_tar]]
                res[oRow][rar_vol] = reactVol[oRow]
                res[oRow][rar_dNTPs] = dicLU_target_dNTPs[res[oRow][rar_tar]]
                res[oRow][rar_dyeConc] = dicLU_target_dyeConc[res[oRow][rar_tar]]
                res[oRow][rar_primer_len_for] = dicLU_target_fw_len[res[oRow][rar_tar]]
                res[oRow][rar_primer_conc_for] = dicLU_target_fw_conc[res[oRow][rar_tar]]
                res[oRow][rar_primer_len_rev] = dicLU_target_rv_len[res[oRow][rar_tar]]
                res[oRow][rar_primer_conc_rev] = dicLU_target_rv_conc[res[oRow][rar_tar]]
                res[oRow][rar_probe1_len] = dicLU_target_probe1_len[res[oRow][rar_tar]]
                res[oRow][rar_probe1_conc] = dicLU_target_probe1_conc[res[oRow][rar_tar]]
                res[oRow][rar_probe2_len] = dicLU_target_probe2_len[res[oRow][rar_tar]]
                res[oRow][rar_probe2_conc] = dicLU_target_probe2_conc[res[oRow][rar_tar]]
                res[oRow][rar_amplicon_len] = dicLU_target_amp_len[res[oRow][rar_tar]]

        if saveRaw:
            rawTable = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl]]]
            for oCol in range(0, spFl[1]):
                rawTable[0].append(oCol + 1)
            for oRow in range(0, spFl[0]):
                rawTable.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar], res[oRow][rar_excl]])
                for oCol in range(0, spFl[1]):
                    rawTable[oRow + 1].append(float(rawFluor[oRow, oCol]))
            finalData["rawData"] = rawTable

        # Count the targets and create the target variables
        # Position 0 is for the general over all window without targets
        vecTarget = np.zeros(spFl[0], dtype=np.int64)
        vecTarget[vecTarget <= 0] = -1
        targetsCount = 1
        tarWinLookup = {}
        tarReverseLookup = {}
        usedTargets = {}
        for oRow in range(0, spFl[0]):
            usedTargets[res[oRow][rar_tar]] = 1.0
            if res[oRow][rar_tar] not in tarWinLookup:
                tarWinLookup[res[oRow][rar_tar]] = targetsCount
                tarReverseLookup[targetsCount] = res[oRow][rar_tar]
                targetsCount += 1
            vecTarget[oRow] = tarWinLookup[res[oRow][rar_tar]]
        upWin = np.zeros(targetsCount, dtype=np.float64)
        lowWin = np.zeros(targetsCount, dtype=np.float64)
        threshold = np.ones(targetsCount, dtype=np.float64)

        # Calculate a report of the supplied data and the limiting factor
        report = ""
        target_limit = {}
        for lu_dye in luDyes:
            report = "The dye \"" + lu_dye.attrib['id'] +"\" is a " + dicLU_dyes[lu_dye.attrib['id']] + " with "
            report += str(dicLU_dye_conc[lu_dye.attrib['id']]) + " &micro;M (&micro;mol/l) dye and "
            report += str(dicLU_dye_dNTPs[lu_dye.attrib['id']]) + " &micro;M (&micro;mol/l) dNTPs.\n"
        for tarID in usedTargets:
            report += "The target \"" + tarID +"\" has a " + str(int(dicLU_target_fw_len[tarID]))
            report += "bp forward primer in " + str(dicLU_target_fw_conc[tarID]) + " nM (nmol/l or "
            report += str(dicLU_target_fw_conc[tarID] * averageVolume / 1000.0) + " pmol in "
            report += str(averageVolume) + " &micro;l), a " + str(int(dicLU_target_rv_len[tarID]))
            report += "bp reverse primer in " + str(dicLU_target_rv_conc[tarID]) + " nM (nmol/l or "
            report += str(dicLU_target_rv_conc[tarID] * averageVolume / 1000.0) + " pmol in "
            report += str(averageVolume) + " &micro;l) resulting in an amplicon of "
            report += str(int(dicLU_target_amp_len[tarID])) + "bp. The exponential phase starts to be limited by "
            # Calculate the limits for each factor
            lowerPrimerConc = dicLU_target_fw_conc[tarID]
            if dicLU_target_rv_conc[tarID] < lowerPrimerConc:
                lowerPrimerConc = dicLU_target_rv_conc[tarID]
            limit_oligo = 0.02 * lowerPrimerConc * 0.000000001 * AVOGADRO * 0.000001
            react_limit_string = "primer concentration "
            target_limit[tarID] = limit_oligo
            if dicLU_targets[tarID] == "non-saturating DNA binding dye":
                limit_dye = 0.1 * 10.4 * dicLU_target_dyeConc[tarID] * 0.000001 * AVOGADRO * 0.000001 / dicLU_target_amp_len[tarID]
                if limit_dye < target_limit[tarID]:
                    react_limit_string = "dye concentration "
                    target_limit[tarID] = limit_dye
            bp_per_amplicon = 2.0 * dicLU_target_amp_len[tarID] - dicLU_target_fw_len[tarID] - dicLU_target_rv_len[tarID]
            limit_dNTPs = 4 * dicLU_target_dNTPs[tarID] * 0.000001 * AVOGADRO * 0.000001 / bp_per_amplicon
            if limit_dNTPs < target_limit[tarID]:
                react_limit_string = "dNTP concentration "
                target_limit[tarID] = limit_dNTPs
            report += react_limit_string + "at " + "{:.3e}".format(int(target_limit[tarID] * averageVolume)) + " copies in "
            report += str(averageVolume) + " &micro;l\n"
        finalData["resultsLinRegPCRReport"] = report

        # Initialization of the error vectors



        vecNoAmplification = np.zeros(spFl[0], dtype=np.bool_)
        vecBaselineError = np.zeros(spFl[0], dtype=np.bool_)
        vecInstableBaseline = np.zeros(spFl[0], dtype=np.bool_)
        vecNoPlateau = np.zeros(spFl[0], dtype=np.bool_)
        vecNoisySample = np.zeros(spFl[0], dtype=np.bool_)
        vecSkipSample = np.zeros(spFl[0], dtype=np.bool_)
        vecShortLogLin = np.zeros(spFl[0], dtype=np.bool_)
        vecCtIsShifting = np.zeros(spFl[0], dtype=np.bool_)
        vecExclEff = np.zeros(spFl[0], dtype=np.bool_)
        vecExclIndivEff = np.zeros(spFl[0], dtype=np.bool_)
        vecTooLowCqEff = np.zeros(spFl[0], dtype=np.bool_)
        vecTooLowCqN0 = np.zeros(spFl[0], dtype=np.bool_)
        vecIsUsedInWoL = np.zeros(spFl[0], dtype=np.bool_)

        # Start and stop cycles of the log lin phase
        rawMinBySD = np.zeros(spFl[0], dtype=np.int64)
        rawMaxBySD = np.zeros(spFl[0], dtype=np.int64)
        maxFD = np.zeros(spFl[0], dtype=np.int64)
        maxSD = np.zeros(spFl[0], dtype=np.int64)
        fluorIncrease = np.zeros(spFl[0], dtype=np.float64)
        rawTD0 = -np.ones(spFl[0], dtype=np.float64)
        check_four_cyc_inc_TD0 = np.zeros(spFl[0], dtype=np.bool_)

        # Baseline correction and PCR eff
        baselineLoopCount = np.zeros(spFl[0], dtype=np.float64)
        baselineError = np.zeros(spFl[0], dtype=np.int64)
        indiv_PCR_Eff = np.ones(spFl[0], dtype=np.float64)


        backgroundNew = -np.ones(spFl[0], dtype=np.float64)


        stopCyc = np.zeros(spFl[0], dtype=np.int64)
        IniStopCyc = np.zeros(spFl[0], dtype=np.int64)
        startCyc = np.zeros(spFl[0], dtype=np.int64)
        startCycFix = np.zeros(spFl[0], dtype=np.int64)

        negShiftBaseline = np.zeros(spFl[0], dtype=np.float64)

        # Initialization of the PCR efficiency vectors
        indiv_thres = np.zeros(spFl[0], dtype=np.float64)

        mean_PCR_Eff = np.zeros(spFl[0], dtype=np.float64)
        mean_indiv_PCR_Eff = np.zeros(spFl[0], dtype=np.float64)
        mean_PCR_Eff_Err = np.zeros(spFl[0], dtype=np.float64)

        iniStart = -np.ones(spFl[0], dtype=np.float64)

        indivCq = -np.ones(spFl[0], dtype=np.float64)
        meanCq = -np.ones(spFl[0], dtype=np.float64)

        td0_mean_thres = np.zeros(spFl[0], dtype=np.float64)
        td0_Cq = -np.ones(spFl[0], dtype=np.float64)
        td0_indiCq = -np.ones(spFl[0], dtype=np.float64)
        td0_meanCq = -np.ones(spFl[0], dtype=np.float64)
        td0_fluor = -np.ones(spFl[0], dtype=np.float64)

        nNulls = np.ones(spFl[0], dtype=np.float64)
        nInclu = np.zeros(spFl[0], dtype=np.int64)
        correl = np.zeros(spFl[0], dtype=np.float64)

        indMeanX = np.zeros(spFl[0], dtype=np.float64)
        indMeanY = np.zeros(spFl[0], dtype=np.float64)

        # Delete
        nNull = np.zeros(spFl[0], dtype=np.float64)
        delNcopy = np.zeros(spFl[0], dtype=np.float64)
        nNull[:] = -1.0
        delNcopy[:] = np.nan

        # Set all to nan
        indMeanX[:] = np.nan
        indMeanY[:] = np.nan


        # Basic Variables
        pointsInWoL = 4

        #####################
        # Fill Matrix Gaps  #
        #####################
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            absMinFluor = np.nanmin(rawFluor)
            if np.isnan(absMinFluor):
                raise RdmlError('Error: Fluorescence data have no valid values.')
            validFloat = np.isfinite(rawFluor)
            gapInWells = {}
            for oRow in range(0, spFl[0]):
                # If a row has no data, use the absolute minimum for all
                if not np.any(validFloat[oRow]):
                    rawFluor[oRow, :] = absMinFluor
                if not np.all(validFloat[oRow]):
                    # Fill gaps at start with first valid value
                    if not validFloat[oRow][0]:
                        lastInvalid = 0
                        while lastInvalid < len(validFloat[oRow]) - 1:
                            lastInvalid += 1
                            if validFloat[oRow][lastInvalid]:
                                break
                        rawFluor[oRow, :lastInvalid] = rawFluor[oRow, lastInvalid]
                        gapInWells[res[oRow][1]] = 0
                    # Fill gaps at end with last valid value
                    if not validFloat[oRow][len(validFloat[oRow]) - 1]:
                        lastInvalid = len(validFloat[oRow])
                        while lastInvalid > 0:
                            lastInvalid -= 1
                            if validFloat[oRow][lastInvalid]:
                                break
                        rawFluor[oRow, (lastInvalid + 1):] = rawFluor[oRow, lastInvalid]
                        gapInWells[res[oRow][1]] = 0
            # Now only internal gaps are remaining
            validFloat = np.isfinite(rawFluor)
            for oRow in range(0, spFl[0]):
                # If a row has no data, use the absolute minimum for all
                if not np.any(validFloat[oRow]):
                    rawFluor[oRow, :] = absMinFluor
                if not np.all(validFloat[oRow]):
                    lastValid = 0
                    gapLen = 0
                    for pos in range(0,len(validFloat[oRow])):
                        if validFloat[oRow][pos]:
                            if gapLen > 0:
                                startFluor = rawFluor[oRow, lastValid]
                                endFluor = rawFluor[oRow, pos]
                                stepFluor = (endFluor - startFluor) / (gapLen + 1)
                                for fillGap in range(1, gapLen + 1):
                                    rawFluor[oRow, lastValid + fillGap] = startFluor + fillGap * stepFluor
                                gapLen = 0
                            else:
                                lastValid = pos
                        else:
                            gapLen += 1
                    gapInWells[res[oRow][1]] = 0
            # This should never happen
            validFloat = np.isfinite(rawFluor)
            if not np.any(validFloat):
                raise RdmlError('Error: Unable to fix all gaps in fluorescence data.')
            orderGapWells = sorted(list(gapInWells.keys()))
            if len(orderGapWells) > 0:
                gapError = "Warning: Extrapolated missing data in well "
                for pos in range(0, len(orderGapWells)):
                    gapError += orderGapWells[pos]
                    if pos < len(orderGapWells) - 1:
                        gapError += ", "
                err = _add_err(err, ";", gapError)

        ################################################
        # Calculate first, second and third derivative #
        ################################################
        nanColAdd =  np.zeros((spFl[0], 1), dtype=np.float64)
        nanColAdd[nanColAdd <= 0.00000001] = np.nan
        tmp = rawFluor - np.roll(rawFluor, 1, axis=1)  # Shift to right
        # tmp = np.roll(np.log10(baselineCorrectedData), 1, axis=1) - np.log10(baselineCorrectedData)  # Shift to right
        rawFlourFD = tmp[:, 1:] # Cq is +0.5
        tmp = rawFlourFD - np.roll(rawFlourFD, 1, axis=1)  # Shift to right
        tmp2 = np.append(nanColAdd, tmp[:, 1:], axis=1)
        rawFlourSD = np.append(tmp2, nanColAdd, axis=1)
        # SD is used as fluorescence mean of 3 cycles
        meanFlourSD = np.append(np.append(nanColAdd, (rawFlourSD[:, 0:-2] + rawFlourSD[:, 1:-1] + rawFlourSD[:, 2:]) / 3.0, axis=1), nanColAdd, axis=1)
        # TD is based on the mean SD
        tmp = meanFlourSD - np.roll(meanFlourSD, 1, axis=1)  # Shift to right
        rawFlourTD = tmp[:, 1:] # Cq is +0.5

        # Min and max fluorescence
        maxBase = np.nanmin(rawFluor, axis=1)
        absMinFluor = np.nanmin(rawFluor)
        absMaxFluor = np.nanmax(rawFluor)
        posFluor = rawFluor + ((absMaxFluor - maxBase) / 100.0 - maxBase)[:, np.newaxis]

        ##################################################
        # Calculate everything what is based on raw data #
        ##################################################
        vecMinFluor = np.nanmin(rawFluor, axis=1)
        # maxFD is FD max + 2 if the fluorescence still increases
        maxFD = np.nanargmax(rawFlourFD, axis=1) + 1  # Cycles so +1 to array
        for row in range(0, spFl[0]):
            if maxFD[row] <= rawFluor.shape[1] - 2:
                # Only add two cycles if there is an increase
                if rawFluor[row, maxFD[row] + 1] > rawFluor[row, maxFD[row]] > rawFluor[row, maxFD[row] - 1]:
                    maxFD[row] += 2
            else:
                maxFD[row] = rawFluor.shape[1]
            maxMeanSD = 0.0
            # maxSD must be before maxFD
            maxSD[row] = rawFluor.shape[1]
            for cyc in range(3, maxFD[row]):
                if not np.isnan(meanFlourSD[row, cyc - 1]) and (meanFlourSD[row, cyc - 1] - maxMeanSD) > 0.0:
                    maxMeanSD = meanFlourSD[row, cyc - 1]
                    maxSD[row] = cyc
            if maxSD[row] + 2 >= rawFluor.shape[1]:
                maxSD[row] = rawFluor.shape[1]
            # rawMinBySD is the Cq to the left with fluorescence decresing
            rawMinBySD[row] = maxSD[row]
            while (rawMinBySD[row] > 3 and rawFluor[row, rawMinBySD[row] - 2] < rawFluor[row, rawMinBySD[row] - 1]):
                rawMinBySD[row] -= 1
            # rawMaxBySD is one Cq to the right if the fluorescence is ascending
            rawMaxBySD[row] = maxSD[row]
            if (rawMaxBySD[row] < rawFluor.shape[1] and rawFluor[row, rawMaxBySD[row] - 1] < rawFluor[row, rawMaxBySD[row]]):
                rawMaxBySD[row] += 1
            # rawTD0 is calculated
            tempStop = min(rawFluor.shape[1] - 1, maxSD[row] + 2)  # Cycles so +1 to array
            tempStop = min(tempStop, rawMaxBySD[row] - 1)  # Cycles so +1 to array
            for cyc in reversed(range(rawMinBySD[row], tempStop)):
                if (rawFlourTD[row, cyc] >= 0.0) and (rawFlourTD[row, cyc + 1] < 0.0):
                    rawTD0[row] = float(cyc) + 1.5 + rawFlourTD[row, cyc] / (rawFlourTD[row, cyc] - rawFlourTD[row, cyc + 1])
                    break
            # Call rawAmplification
            if int(rawTD0[row]) - rawMinBySD[row] + 1 > 3:
                check_four_cyc_inc_TD0[row] = True
            else:
                vecSkipSample[row] = True
        #        rawTD0[row] = -1.0

        ########################################
        # Correct negative fluorescence values #
        ########################################
        posFluor = rawFluor.copy()
        # There should be no negative values in uncorrected raw data
        absMinFluor = np.nanmin(rawFluor)
        absMaxFluor = np.nanmax(rawFluor)
        if absMinFluor < 0.00000001:
            err = _add_err(err, ";", "Error: Fluorescence data have negative values. Use raw data without baseline correction! ")
            err = _add_err(err, "", "Baseline corrected data not using a constant factor will result in wrong PCR efficiencies!")
            vecMinFluor = np.nanmin(rawFluor, axis=1)
            for row in range(0, spFl[0]):
                if vecMinFluor[row] < 0.00000001:
                    negShiftBaseline[row] = (absMaxFluor - vecMinFluor[row]) / 100.0 - vecMinFluor[row]
                    for col in range(0, spFl[1]):
                        posFluor[row, col] += negShiftBaseline[row]

        ######################################################
        # Baseline correction and PCR efficiency calculation #
        ######################################################
        # Set a background for failing samples
        vecBackground = np.nanmin(posFluor, axis=1) * 0.99
        for row in range(0, spFl[0]):
            # start and maxSD in Cycles, so +1 to array.
            # Initialize start values
            curBase = 0.0
            baselineError[row] = 5
            stop = int(rawTD0[row])
            # At minimum 4 cycles are required
            if stop < 4:
                baselineError[row] = 1
                continue
            # At maximum 9 cycles should be inluded
            start = max(rawMinBySD[row] + 1, stop - 9 + 1)
            # Now check that 4 cycles are left
            if stop - start + 1 < 4:
                baselineError[row] = 2
                continue
            selSection = posFluor[row][start - 1:stop]
            # Check if the baseline is within the range
            resInc = _lrp_testSlopes2(selSection)
            if resInc[0] > resInc[1] :
                baselineError[row] = 3
                continue
            testSelection = selSection - posFluor[row][start - 1] + 0.00000001
            resInc = _lrp_testSlopes2(testSelection)
            if resInc[0] < resInc[1] :
                baselineError[row] = 4
                continue
            # Half the area and continue search in the good part
            curBase = posFluor[row][start - 1] * 0.5
            step = curBase * 0.5
            currentSection = selSection - curBase
            for loop in range (0, 100):
                resInc = _lrp_testSlopes2(currentSection)
                if (np.abs(resInc[0] - resInc[1]) < 0.00001):
                    vecBackground[row] = curBase
                    indiv_PCR_Eff[row] = np.power(10, resInc[1])
                    baselineLoopCount[row] = loop
                    baselineError[row] = 0
                    break
                if resInc[0] > resInc[1]:
                    curBase -= step
                else:
                    curBase += step
                step *= 0.5
                currentSection = selSection - curBase

        # Do the baseline correction
        baselineCorrectedData = posFluor - vecBackground[:, np.newaxis]

        # Calculate the baseline corrected fluorescence for the TD0
        for row in range(0, spFl[0]):
            if check_four_cyc_inc_TD0[row]:
                low = int(np.floor(rawTD0[row])) - 1
                high = int(np.ceil(rawTD0[row])) - 1
                td0_fluor[row] = np.power(10, np.log10(baselineCorrectedData[row, low]) + (np.log10(baselineCorrectedData[row, high]) - np.log10(baselineCorrectedData[row, low])) * (rawTD0[row] - np.floor(rawTD0[row] - 1)))

        # Calulate plateau baseline ratio
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            meanPlateau = np.nanmean(baselineCorrectedData[:, -5:], axis=1)

        if saveBaslineCorr:
            rawTable = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl]]]
            for oCol in range(0, spFl[1]):
                rawTable[0].append(oCol + 1)
            for oRow in range(0, spFl[0]):
                rawTable.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar], res[oRow][rar_excl]])
                for oCol in range(0, spFl[1]):
                    rawTable[oRow + 1].append(float(baselineCorrectedData[oRow, oCol]))
            finalData["baselineCorrectedData"] = rawTable


        # Negative controls should not be part of the mean calculations
        for row in range(0, spFl[0]):
            if res[row][rar_sample_type] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                vecSkipSample[row] = True


        # TODO: Remove
        backgroundNew = vecBackground.copy()
        WoL_PCR_Eff = indiv_PCR_Eff.copy()


        




        #######################
        # Baseline correction #
        #######################
        start_time = datetime.datetime.now()

        ###########################################################################
        # First quality check : Is there enough amplification during the reaction #
        ###########################################################################

        ffocuss = "xxx"

        # Slope calculation per react/target - the intercept is never used for now
        baseCorFluor = posFluor.copy()
        rawMod = posFluor.copy()
        # print(rawFluor)

        slopeAmp = _lrp_linReg(adp_cyc_min, np.log10(rawMod))

        # Calculate the minimum of fluorescence values per react/target, store it as background
        # and substract it from the raw fluorescence values
        vecMinFluor = np.nanmin(rawMod, axis=1)
        vecBackground = 0.99 * vecMinFluor
        vecDefBackgrd = vecBackground.copy()
        minCorFluor = rawMod - vecBackground[:, np.newaxis]
        # print(vecMinFluor)
        # print(minCorFluor)

        minFluCount = np.ones(minCorFluor.shape, dtype=np.int64)
        minFluCountSum = np.sum(minFluCount, axis=1)
        minSlopeAmp = _lrp_linReg(adp_cyc_min, np.log10(minCorFluor))
        
        for oRow in range(0, spFl[0]):
            # Check to detect the negative slopes and the PCR reactions that have an
            # amplification less than seven the minimum fluorescence
            if slopeAmp[oRow] < 0 or minSlopeAmp[oRow] < (np.log10(7.0) / minFluCountSum[oRow]):
                vecNoAmplification[oRow] = True

            # There must be an increase in fluorescence after the amplification.
            # TODO take the min of all
            if minCorFluor[oRow, -1] / np.nanmean(minCorFluor[oRow, 0:10]) < 7:
                vecNoAmplification[oRow] = True

            if not vecNoAmplification[oRow]:
                IniStopCyc[oRow] = stopCyc[oRow] = maxSD[oRow]
                [startCyc[oRow], startCycFix[oRow]] = _lrp_findStartCyc(minCorFluor, oRow, stopCyc[oRow])
                [iniStart[oRow], blaaa] = _lrp_findStartCyc(rawMod, oRow, stopCyc[oRow])
            else:
                vecSkipSample[oRow] = True
                stopCyc[oRow] = minCorFluor.shape[1]
                IniStopCyc[oRow] = stopCyc[oRow]
                startCyc[oRow] = 1
                startCycFix[oRow] = 1

            if 0:
                print(maxSD[oRow])
                print(rawMaxBySD[oRow])
                print(stopCyc[oRow])
                print("        ")

            # iniStart[oRow] = startCyc[oRow]

            # Get the positions ignoring nan values
            posCount = 0
            posMinOne = 0
            posMinTwo = 0
            for realPos in range(stopCyc[oRow] - 2, 0, -1):
                if not np.isnan(minCorFluor[oRow, realPos - 1]):
                    if posCount == 0:
                        posMinOne = realPos + 1
                    if posCount > 0:
                        posMinTwo = realPos + 1
                        break
                    posCount += 1

            if not (minCorFluor[oRow, stopCyc[oRow] - 1] > minCorFluor[oRow, posMinOne - 1] > minCorFluor[oRow, posMinTwo - 1]):
                vecNoAmplification[oRow] = True
                vecSkipSample[oRow] = True

            if vecNoAmplification[oRow] or vecBaselineError[oRow] or stopCyc[oRow] == minCorFluor.shape[1]:
                vecNoPlateau[oRow] = True

        # Set an initial window already for WOL calculation
        lastCycMeanMax = _lrp_lastCycMeanMax(minCorFluor, vecSkipSample, vecNoPlateau)
        upWin[0] = 0.1 * lastCycMeanMax
        lowWin[0] = 0.1 * lastCycMeanMax / 16.0

        

        # Write over
        vecBackground = backgroundNew
        baseCorFluor = baselineCorrectedData
        meanTarEff = {}


        vecSkipSample[vecExcludedByUser] = True
        # Update the window
        lastCycMeanMax = _lrp_lastCycMeanMax(baseCorFluor, vecSkipSample, vecNoPlateau)
        upWin[0] = 0.1 * lastCycMeanMax
        lowWin[0] = 0.1 * lastCycMeanMax / 16.0
        maxFluorTotal = np.nanmax(baseCorFluor)
        minFluorTotal = np.nanmin(baseCorFluor)
        if minFluorTotal < maxFluorTotal / 10000:
            minFluorTotal = maxFluorTotal / 10000

        # Fixme: Per group
        # CheckNoisiness
        skipGroup = False
        maxLim = _lrp_meanStopFluor(baseCorFluor, None, None, stopCyc, vecSkipSample, vecNoPlateau)
        if maxLim > 0.0:
            maxLim = np.log10(maxLim)
        else:
            skipGroup = True
        checkMeanEff = 1.0

        if not skipGroup:
            step = pointsInWoL * _lrp_logStepStop(baseCorFluor, None, [], stopCyc, vecSkipSample, vecNoPlateau)
            upWin, lowWin = _lrp_setLogWin(None, maxLim, step, upWin, lowWin, maxFluorTotal, minFluorTotal)
            # checkBaseline = np.log10(0.5 * np.round(1000 * np.power(10, upWin[0])) / 1000)
            _unused, _unused2, tempPcrEff, _unused3, _unused4, _unused5 = _lrp_allParamInWindow(baseCorFluor,
                                                                                                None, [],
                                                                                                indMeanX, indMeanY,
                                                                                                WoL_PCR_Eff, nNulls,
                                                                                                nInclu, correl,
                                                                                                startCyc, startCycFix, stopCyc, 
                                                                                                upWin, lowWin,
                                                                                                vecNoAmplification,
                                                                                                vecBaselineError)
            checkMeanEff, _unused = _lrp_meanPcrEff(None, [], tempPcrEff, vecSkipSample, vecNoPlateau, vecShortLogLin)
            if checkMeanEff < 1.001:
                skipGroup = True

        if not skipGroup:
            foldWidth = np.log10(np.power(checkMeanEff, pointsInWoL))
            upWin, lowWin = _lrp_setLogWin(None, maxLim, foldWidth, upWin, lowWin, maxFluorTotal, minFluorTotal)
            # compare to Log(1.01*lowLim) to compensate for
            # the truncation in cuplimedit with + 0.0043
            lowLim = maxLim - foldWidth + 0.0043
            for oRow in range(0, spFl[0]):
                if not vecSkipSample[oRow]:
                    startWinCyc, stopWinCyc, _unused = _lrp_startStopInWindow(baseCorFluor, oRow, startCyc, startCycFix, stopCyc, upWin[0], lowWin[0])
                    minStartCyc = startWinCyc - 1
                    # Handle possible NaN
                    while np.isnan(baseCorFluor[oRow, minStartCyc - 1]) and minStartCyc > 1:
                        minStartCyc -= 1
                    minStopCyc = stopWinCyc - 1
                    while np.isnan(baseCorFluor[oRow, minStopCyc - 1]) and minStopCyc > 2:
                        minStopCyc -= 1

                    minStartFlour = baseCorFluor[oRow, minStartCyc - 1]
                    if np.isnan(minStartFlour):
                        minStartFlour = 0.00001

                    startStep = np.log10(baseCorFluor[oRow, startWinCyc - 1]) - np.log10(minStartFlour)
                    stopStep = np.log10(baseCorFluor[oRow, stopWinCyc - 1]) - np.log10(baseCorFluor[oRow, minStopCyc - 1])
                    if (np.log10(minStartFlour) > lowLim and not
                            ((minStartFlour < baseCorFluor[oRow, startWinCyc - 1] and startStep < 1.2 * stopStep) or
                             (startWinCyc - minStartCyc > 1.2))):
                        vecNoisySample[oRow] = True
                        vecSkipSample[oRow] = True



        if timeRun:
            stop_time = datetime.datetime.now() - start_time
            print("Done Baseline: " + str(stop_time) + "sec")

        ###########################################################
        # Calculation of the Window of Linearity (WOL) per target #
        ###########################################################

        # Set a starting window for all groups
        for tar in range(1, targetsCount):
            upWin[tar] = upWin[0]
            lowWin[tar] = lowWin[0]

        for oRow in range(0, spFl[0]):
            if vecNoAmplification[oRow] or vecBaselineError[oRow] or stopCyc[oRow] == spFl[1]:
                vecNoPlateau[oRow] = True
            else:
                vecNoPlateau[oRow] = False

        for tar in range(1, targetsCount):
            indMeanX, indMeanY, WoL_PCR_Eff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL = _lrp_setWoL(baseCorFluor, tar, vecTarget, pointsInWoL,
                                                                                                                       indMeanX, indMeanY, WoL_PCR_Eff, nNulls, nInclu,
                                                                                                                       correl, upWin, lowWin, maxFluorTotal,
                                                                                                                       minFluorTotal, stopCyc, startCyc, startCycFix, threshold,
                                                                                                                       vecNoAmplification, vecBaselineError,
                                                                                                                       vecSkipSample, vecNoPlateau, vecShortLogLin,
                                                                                                                       vecIsUsedInWoL)
            indMeanX, indMeanY, WoL_PCR_Eff, nNulls, nInclu, correl, upWin, lowWin, threshold, vecIsUsedInWoL, vecNoPlateau = _lrp_assignNoPlateau(baseCorFluor, tar, vecTarget,
                                                                                                                                              pointsInWoL, indMeanX, indMeanY,
                                                                                                                                              WoL_PCR_Eff, nNulls, nInclu, correl,
                                                                                                                                              upWin, lowWin, maxFluorTotal,
                                                                                                                                              minFluorTotal, stopCyc, startCyc, startCycFix,
                                                                                                                                              threshold, vecNoAmplification,
                                                                                                                                              vecBaselineError, vecSkipSample,
                                                                                                                                              vecNoPlateau, vecShortLogLin,
                                                                                                                                              vecIsUsedInWoL)

        for oRow in range(0, spFl[0]):
            stopTmp = stopCyc[oRow] + 2
            maxLen = len(rawFlourTD[oRow])
            if maxLen < stopTmp + 1:
                stopTmp = maxLen - 1
            for cyc in reversed(range(startCyc[oRow], stopTmp)):
              #  if oRow in [242, 277]:
              #      print(str(stopTmp) + " - " + str(cyc) + " - "  + str(rawFlourTD[oRow, cyc]) + " - " + str(rawFlourTD[oRow, cyc + 1]) + " - " )
                if np.isnan(np.isnan(rawFlourTD[oRow, cyc])):
                    continue
                if np.isnan(np.isnan(rawFlourTD[oRow, cyc + 1])):
                    continue
                if (rawFlourTD[oRow, cyc] >= 0.0) and (rawFlourTD[oRow, cyc + 1] < 0.0):
                    td0_Cq[oRow] = rawTD0[oRow]
                    low = int(np.floor(td0_Cq[oRow])) - 1
                    high = int(np.ceil(td0_Cq[oRow])) - 1
                    td0_fluor[oRow] = np.power(10, np.log10(baselineCorrectedData[oRow, low]) + ( np.log10(baselineCorrectedData[oRow, high]) - np.log10(baselineCorrectedData[oRow, low])) * (td0_Cq[oRow] - np.floor(td0_Cq[oRow])))
                    break
    #        print(res[oRow][rar_well] + " Inc: " + str(fluorIncrease[oRow]) + " Start: " + str(rawMinBySD[oRow]) + " Stop: " + str(maxSD[oRow]) + " TD0: " + str(rawTD0[oRow]) + " - " + str(td0_Cq[oRow]))

        # Without window of lineartity
       # indiv_PCR_Eff = saveEff



          #      for bla in dddddd[oRow]:
           #         pTab[oRow + 1].append(str(bla))
                
            if 0:
                ww = open("/home/untergasser/code/rdml/tools/server/rdmlpython/test/temp_td0.csv", "w")
                for row in range(0, len(pTab)):
                    for col in range(0, len(pTab[row])):
                        ww.write(str(pTab[row][col]) + "\t")
                    ww.write("\n")
                ww.close()

                for tar in eff_perTarget:
                    print(tar)
                    for pos in range(0, len(eff_perTarget[tar]["raw"])):
                        print("  " + str(eff_perTarget[tar]["raw"][pos]) + "  " + str(eff_wol_perTarget[tar]["raw"][pos]) + "  " + str(corel_perTarget[tar]["raw"][pos]))
                    eff_perTarget[tar]["mean"] = np.mean(eff_perTarget[tar]["raw"])
                    eff_perTarget[tar]["std"] = np.std(eff_perTarget[tar]["raw"], ddof=1)
                    eff_wol_perTarget[tar]["mean"] = np.mean(eff_wol_perTarget[tar]["raw"])
                    eff_wol_perTarget[tar]["std"] = np.std(eff_wol_perTarget[tar]["raw"], ddof=1)
                for tar in saveTD0:
                    saveTD0[tar]["mean"] = np.mean(saveTD0[tar]["raw"])
                    saveTD0[tar]["std"] = np.std(saveTD0[tar]["raw"], ddof=1)
                    saveTD0fluor[tar]["mean"] = np.mean(saveTD0fluor[tar]["raw"])
                    saveTD0fluor[tar]["std"] = np.std(saveTD0fluor[tar]["raw"], ddof=1)

            if 0:
                for tar in saveTD0:
                    if tar in saveTD0:
                        print(tar + " " + "{:.4f}".format(saveTD0[tar]["mean"]) + " (" + "{:.4f}".format(saveTD0[tar]["std"]) + ") " + "{:.4f}".format(saveTD0fluor[tar]["mean"]) + " (" + "{:.4f}".format(saveTD0fluor[tar]["std"]) + ")")

            if 0:
                for tar in eff_orddr:
                    dil = -1.0
                    if tar in effDiutions:
                        dil = effDiutions[tar]
                    if tar in eff_perTarget:
                        print(tar + " " + "{:.4f}".format(dil) + " " + "{:.4f}".format(dil - eff_perTarget[tar]["mean"]) + "   " + "{:.4f}".format(eff_perTarget[tar]["mean"]) + " (" + "{:.4f}".format(eff_perTarget[tar]["std"]) + ")" 
                                + "   " + "{:.4f}".format(eff_wol_perTarget[tar]["mean"]) + " (" + "{:.4f}".format(eff_wol_perTarget[tar]["std"]) + ")")


        # Use WOL or not 
        # WoL_PCR_Eff = indiv_PCR_Eff


        # Calc Geomean fluor TD0:
        for tar in range(1, targetsCount):
            fff_num = 0
            fff_sum = 0.0
            fff_res = -1.0
            for oRow in range(0, spFl[0]):
                if res[oRow][rar_sample_type] not in ["std", "pos", "unkn"]:
                    continue
                if vecTarget[oRow] == tar:
                    if td0_fluor[oRow] > 0.0:
                        fff_num += 1
                        fff_sum += td0_fluor[oRow]  # math.log(td0_fluor[oRow])
            if fff_num > 0:
                fff_res = fff_sum / fff_num  # math.exp(fff_sum / fff_num)
            for oRow in range(0, spFl[0]):
                if vecTarget[oRow] == tar:
                    td0_mean_thres[oRow] = fff_res

        # TODO: Fix this
        for oRow in range(0, spFl[0]):
            if td0_mean_thres[oRow] > 0.0:
                td0_indiCq[oRow] = indMeanX[oRow] + (np.log10(td0_mean_thres[oRow]) - indMeanY[oRow]) / np.log10(WoL_PCR_Eff[oRow])  # td0_indiCq[oRow]  td0_Cq[oRow]

        # TODO Maybe an error
    #   for oRow in range(0, spFl[0]):
      #      if not vecNoAmplification[oRow] and td0_indiCq[oRow] < 0.0:
      #          print("TD0 not found: " + str(oRow) + " Well: " + res[oRow][rar_well] + " Base: " + str(vecBaselineError[oRow]))

        # Median values calculation
        vecSkipSample_Plat = vecSkipSample.copy()
        vecSkipSample_Plat[vecNoPlateau] = True

        # Common threshold calculation
        comThreshSum = 0.0
        comThreshNum = 0
        for oRow in range(0, spFl[0]):
            if not vecSkipSample[oRow]:
                if not vecNoAmplification[oRow]:
                    if res[oRow][rar_sample_type] in ["std", "pos", "unkn"]:
                        comThreshSum += math.log(threshold[vecTarget[oRow]])
                        comThreshNum += 1
        if comThreshNum > 0:
            threshold[0] = math.exp(comThreshSum / comThreshNum)
        else:
            logThreshold = np.log10(threshold[1:])
            threshold[0] = np.power(10, np.mean(logThreshold))

        # Create the warnings for the different chemistries
        # Chem Arr     0     1     2     3     4     5     6     7     8     9    10
        critCqEff = [28.0, 28.0, 19.0, 16.0, 14.0, 12.0, 11.0, 11.0, 10.0, 10.0,  9.0]  # For error Eff < 0.01
        critCqN0 = [40.0, 40.0, 27.0, 19.0, 16.0, 13.0, 12.0, 11.0, 10.0,  9.0,  9.0]  # For bias N0 < 0.95
        for oRow in range(0, spFl[0]):
            if res[oRow][rar_tar_chemistry] in ["hydrolysis probe", "labelled reverse primer", "DNA-zyme probe"]:
                critCqOffset = 0.0
                if (res[oRow][rar_tar_chemistry] == "labelled reverse primer" and
                        res[oRow][rar_sample_nucleotide] == "ds"):
                    critCqOffset = 1.0
                if (res[oRow][rar_tar_chemistry] == "DNA-zyme probe" and
                        res[oRow][rar_sample_nucleotide] == "ds"):
                    critCqOffset = 4.0
                if (res[oRow][rar_tar_chemistry] == "DNA-zyme probe" and
                        res[oRow][rar_sample_nucleotide] == "ss"):
                    critCqOffset = 6.0
                if (not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0001 and
                        threshold[vecTarget[oRow]] > 0.0001 and not (vecNoAmplification[oRow] or vecBaselineError[oRow])):
                    effIndex = int(np.trunc(10 * WoL_PCR_Eff[oRow] + 1 - 10))
                    if effIndex < 0:
                        effIndex = 0
                    if effIndex > 10:
                        effIndex = 10
                    tempCq_Grp = indMeanX[oRow] + (np.log10(threshold[0]) - indMeanY[oRow]) / np.log10(WoL_PCR_Eff[oRow])
                    if tempCq_Grp > 0.0:
                        if tempCq_Grp < (critCqEff[effIndex] + critCqOffset):
                            vecTooLowCqEff[oRow] = True
                        if tempCq_Grp < (critCqN0[effIndex] + critCqOffset):
                            vecTooLowCqN0[oRow] = True

        pcreff_NoNaN = WoL_PCR_Eff.copy()
        pcreff_NoNaN[np.isnan(WoL_PCR_Eff)] = 0.0
        for tar in range(1, targetsCount):
            # Calculating all choices takes less time then to recalculate
            pcreff_Skip = WoL_PCR_Eff.copy()
            pcreff_Skip[vecTooLowCqEff] = np.nan
            pcreff_Skip[vecSkipSample] = np.nan
            pcreff_Skip[pcreff_NoNaN < 1.001] = np.nan
            pcreff_Skip[~(vecTarget == tar)] = np.nan

            pcreff_indiv_Skip = indiv_PCR_Eff.copy()
            pcreff_indiv_Skip[vecTooLowCqEff] = np.nan
            pcreff_indiv_Skip[vecSkipSample] = np.nan
            pcreff_indiv_Skip[pcreff_NoNaN < 1.001] = np.nan
            pcreff_indiv_Skip[~(vecTarget == tar)] = np.nan

            if excludeNoPlateau is True:
                pcreff_Skip[vecNoPlateau] = np.nan

            if excludeEfficiency == "mean":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pcreffMedian_Skip = np.nanmedian(pcreff_Skip)
                for oRow in range(0, spFl[0]):
                    if tar == vecTarget[oRow]:
                        if not np.isnan(WoL_PCR_Eff[oRow]):
                            if (np.isnan(pcreffMedian_Skip) or
                                    not (pcreffMedian_Skip - pcrEfficiencyExl <= WoL_PCR_Eff[oRow] <= pcreffMedian_Skip + pcrEfficiencyExl)):
                                vecExclEff[oRow] = True
                pcreff_Skip[vecExclEff] = np.nan

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pcreffMedian_Skip = np.nanmedian(pcreff_Skip)
                for oRow in range(0, spFl[0]):
                    if tar is None or tar == vecTarget[oRow]:
                        if not np.isnan(WoL_PCR_Eff[oRow]):
                            if (np.isnan(pcreffMedian_Skip) or
                                    not (pcreffMedian_Skip - pcrEfficiencyExl <= WoL_PCR_Eff[oRow] <= pcreffMedian_Skip + pcrEfficiencyExl)):
                                vecExclEff[oRow] = True
                            else:
                                vecExclEff[oRow] = False
                        else:
                            vecExclEff[oRow] = True
                pcreff_Skip[vecExclEff] = np.nan

                if 0:
                    with warnings.catch_warnings():
                        print(pcreff_indiv_Skip)
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        pcreffMedian_Skip = np.nanmedian(pcreff_indiv_Skip)
                    for oRow in range(0, spFl[0]):
                        if tar == vecTarget[oRow]:
                            if not np.isnan(indiv_PCR_Eff[oRow]):
                                print(str(pcreffMedian_Skip - pcrEfficiencyExl) + " " + str(indiv_PCR_Eff[oRow]) + " " + str(pcreffMedian_Skip + pcrEfficiencyExl))
                                if (np.isnan(pcreffMedian_Skip) or
                                        not (pcreffMedian_Skip - pcrEfficiencyExl <= indiv_PCR_Eff[oRow] <= pcreffMedian_Skip + pcrEfficiencyExl)):
                                    vecExclIndivEff[oRow] = True
                                    print("Removed")
                    pcreff_indiv_Skip[vecExclIndivEff] = np.nan

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    pcreffMedian_Skip = np.nanmedian(pcreff_indiv_Skip)
                for oRow in range(0, spFl[0]):
                    if tar is None or tar == vecTarget[oRow]:
                        if not np.isnan(indiv_PCR_Eff[oRow]):
                            if (np.isnan(pcreffMedian_Skip) or
                                    not (pcreffMedian_Skip - pcrEfficiencyExl <= indiv_PCR_Eff[oRow] <= pcreffMedian_Skip + pcrEfficiencyExl)):
                                vecExclIndivEff[oRow] = True
                            else:
                                vecExclIndivEff[oRow] = False
                        else:
                            vecExclIndivEff[oRow] = True
                pcreff_indiv_Skip[vecExclIndivEff] = np.nan

            if excludeEfficiency == "outlier":
                vecExclEff[_lrp_removeOutlier(pcreff_Skip, vecNoPlateau)] = True
                vecExclIndivEff[_lrp_removeOutlier(pcreff_indiv_Skip, vecNoPlateau)] = True
                pcreff_Skip[vecExclEff] = np.nan
                pcreff_indiv_Skip[vecExclIndivEff] = np.nan

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                tempMeanEff_Skip = np.nanmean(pcreff_Skip)
                tempMeanIndiEff_Skip = np.nanmean(pcreff_indiv_Skip)
                tempStdEff_Skip = np.nanstd(pcreff_Skip)

            for oRow in range(0, spFl[0]):
                if tar == vecTarget[oRow]:
                    mean_PCR_Eff[oRow] = tempMeanEff_Skip
                    mean_indiv_PCR_Eff[oRow] = tempMeanIndiEff_Skip
                    mean_PCR_Eff_Err[oRow] = tempStdEff_Skip
                    meanTarEff[tarReverseLookup[tar]] = tempMeanEff_Skip

                    # Delete
                    if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0001 and threshold[tar] > 0.0001 and not (vecNoAmplification[oRow] or vecBaselineError[oRow]):
                        indivCq[oRow] = indMeanX[oRow] + (np.log10(threshold[0]) - indMeanY[oRow]) / np.log10(WoL_PCR_Eff[oRow])

                    if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0:
                        if not np.isnan(mean_PCR_Eff[oRow]) and mean_PCR_Eff[oRow] > 1.001:
                            meanCq[oRow] = indMeanX[oRow] + (np.log10(threshold[0]) - indMeanY[oRow]) / np.log10(mean_PCR_Eff[oRow])


                    if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0 and 0.0 < indivCq[oRow] < 2 * spFl[1]:
                        if not np.isnan(mean_PCR_Eff[oRow]) and mean_PCR_Eff[oRow] > 1.001:
                            nNull[oRow] = threshold[0] / np.power(mean_PCR_Eff[oRow], meanCq[oRow])


                    # Correction of the different chemistries
             #       cqCorrection = 0.0
             #       if res[oRow][rar_tar_chemistry] in ["hydrolysis probe", "labelled reverse primer", "DNA-zyme probe"]:
             #           cqCorrection = -1.0

              #      if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0001 and threshold[tar] > 0.0001 and not (vecNoAmplification[oRow] or vecBaselineError[oRow]):
              #          if res[oRow][rar_tar_chemistry] == "DNA-zyme probe":
             #               cqCorrection = -1.0 + np.log10(1 / (1 - (1 / WoL_PCR_Eff[oRow]))) / np.log10(WoL_PCR_Eff[oRow])
              #          indivCq[oRow] = indMeanX[oRow] + (np.log10(threshold[0]) - indMeanY[oRow]) / np.log10(WoL_PCR_Eff[oRow]) + cqCorrection
             #           indivCq_Grp[oRow] = indMeanX[oRow] + (np.log10(threshold[tar]) - indMeanY[oRow]) / np.log10(WoL_PCR_Eff[oRow]) + cqCorrection

             #       if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0:
             #           if not np.isnan(mean_PCR_Eff[oRow]) and mean_PCR_Eff[oRow] > 1.001:
              #              if res[oRow][rar_tar_chemistry] == "DNA-zyme probe":
              #                  cqCorrection = -1.0 + np.log10(1 / (1 - (1 / mean_PCR_Eff[oRow]))) / np.log10(mean_PCR_Eff[oRow])
              #              meanCq_Skip[oRow] = indMeanX[oRow] + (np.log10(threshold[0]) - indMeanY[oRow]) / np.log10(mean_PCR_Eff[oRow]) + cqCorrection

              #      if not np.isnan(WoL_PCR_Eff[oRow]) and WoL_PCR_Eff[oRow] > 1.0 and 0.0 < indivCq[oRow] < 2 * spFl[1]:
              #          if not np.isnan(mean_PCR_Eff[oRow]) and mean_PCR_Eff[oRow] > 1.001:
               #             meanNnull_Skip[oRow] = threshold[0] / np.power(mean_PCR_Eff[oRow], meanCq_Skip[oRow])

        if 0:
            for tar in eff_orddr:
                dil = -1.0
                if tar in effDiutions:
                    dil = effDiutions[tar]
                if tar in eff_perTarget:
                    print(tar + " " + "{:.4f}".format(dil) + "   " + "{:.4f}".format(meanTarEff[tar]) + "    " + "{:.4f}".format(eff_perTarget[tar]["mean"]) + " (" + "{:.4f}".format(eff_perTarget[tar]["std"]) + ")")


        if 0:
            saveTD0 = {}
            saveTD0fluor = {}
            eff_perTarget = {}
            eff_wol_perTarget = {}
            corel_perTarget = {}


            pTab = [["run", header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl], "start", "end", "TD0", "TD0 fluor", "TD0 end", "4cycInc", "BG", "Baseline Err", "Dilution Eff", "indiv PCR Eff", "WoL PCR Eff"]]
            for oRow in range(0, spFl[0]):
                if res[oRow][rar_sample_type] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                    continue
                dilStd = ""
                if res[oRow][rar_tar] in effDiutions:
                    dilStd = str(effDiutions[res[oRow][rar_tar]])
                pTab.append([self["id"] , res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar], res[oRow][rar_excl], rawMinBySD[oRow], maxSD[oRow], rawTD0[oRow], 
                             td0_fluor[oRow], int(rawTD0[oRow]), check_four_cyc_inc_TD0[oRow], backgroundNew[oRow], baselineError[oRow], dilStd, indiv_PCR_Eff[oRow], WoL_PCR_Eff[oRow]])
                if baselineError[oRow] == 0:
                    if res[oRow][rar_tar] not in eff_perTarget:
                        eff_perTarget[res[oRow][rar_tar]] = {}
                        eff_perTarget[res[oRow][rar_tar]]["raw"] = []
                        eff_wol_perTarget[res[oRow][rar_tar]] = {}
                        eff_wol_perTarget[res[oRow][rar_tar]]["raw"] = []
                        corel_perTarget[res[oRow][rar_tar]] = {}
                        corel_perTarget[res[oRow][rar_tar]]["raw"] = []

                    eff_perTarget[res[oRow][rar_tar]]["raw"].append(indiv_PCR_Eff[oRow])
                    eff_wol_perTarget[res[oRow][rar_tar]]["raw"].append(WoL_PCR_Eff[oRow])
                    corel_perTarget[res[oRow][rar_tar]]["raw"].append(correl[oRow] * correl[oRow])

                    if res[oRow][rar_tar] not in saveTD0:
                        saveTD0[res[oRow][rar_tar]] = {}
                        saveTD0[res[oRow][rar_tar]]["raw"] = []
                        saveTD0fluor[res[oRow][rar_tar]] = {}
                        saveTD0fluor[res[oRow][rar_tar]]["raw"] = []
                    saveTD0[res[oRow][rar_tar]]["raw"].append(rawTD0[oRow])
                    saveTD0fluor[res[oRow][rar_tar]]["raw"].append(td0_fluor[oRow])




        effDiutions = {"FSTL_1_*_259": 1.931, "FSTL_1_*_259 SF": 1.9044, "FSTL_1_A_047": 1.9705, "FSTL_1_A_047 SF": 1.9836, "FSTL_1_B_042": 1.956, "FSTL_1_D_105": 1.9327, 
                    "FSTL_1_D_105 SF": 1.9765, "FSTL_1_E_097": 1.9604, "FSTL_1_E_097 SF": 1.9456, "FSTL_1_F_109": 1.9253, "FSTL_1_F_109 LC": 2.1111, "FSTL_1_F_109 SF": 1.9012, 
                    "FSTL_1_H_201": 1.8712, "FSTL_1_H_201 LC": 1.9652, "FSTL_1_H_201 SF": 1.8351, "FSTL_1_I_204": 1.9444, "FSTL_1_I_204 LC": 1.9291, "FSTL_1_I_204 SF": 1.9129, 
                    "FSTL_1_K_219": 1.9418, "FSTL_1_K_219 LC": 1.9682, "FSTL_1_K_219 SF": 1.9061, "FSTL_1_L_412": 1.874, "FSTL_1_L_412 SF": 1.858, "FSTL_1_M_398": 2.2317, 
                    "FSTL_1_M_398 LC": 1.9586, "FSTL_1_M_398 SF": 1.8453, "FSTL_1_O_417 LC": 1.9948, "FSTL_1_O_417 SF": 1.871, "FSTL_1_P_820": 2.0972, "FSTL_1_P_820 LC": 1.8312, 
                    "FSTL_1_P_820 SF": 1.7665, "FSTL_1_T_770 SF": 1.7915}

        eff_orddr = ["FSTL_1_A_047 100nM", "FSTL_1_A_047", "FSTL_1_A_047 750nM", "FSTL_1_B_042 100nM", "FSTL_1_B_042", "FSTL_1_B_042 750nM", "FSTL_1_D_105 100nM", "FSTL_1_D_105", "FSTL_1_D_105 750nM", 
                     "FSTL_1_E_097 100nM", "FSTL_1_E_097", "FSTL_1_E_097 750nM", "FSTL_1_F_109 100nM", "FSTL_1_F_109", "FSTL_1_F_109 750nM", "FSTL_1_H_201 100nM", "FSTL_1_H_201", "FSTL_1_H_201 750nM", 
                     "FSTL_1_I_204 100nM", "FSTL_1_I_204", "FSTL_1_I_204 750nM", "FSTL_1_K_219 100nM", "FSTL_1_K_219", "FSTL_1_K_219 750nM", "FSTL_1_*_259 100nM", "FSTL_1_*_259", "FSTL_1_*_259 750nM", 
                     "FSTL_1_L_412 100nM", "FSTL_1_L_412", "FSTL_1_L_412 750nM", "FSTL_1_M_398 100nM", "FSTL_1_M_398", "FSTL_1_M_398 750nM", "FSTL_1_P_820 100nM", "FSTL_1_P_820", "FSTL_1_P_820 750nM", 
                     "FSTL_1_A_047 100nM SF", "FSTL_1_A_047 SF", "FSTL_1_A_047 750nM SF", "FSTL_1_D_105 100nM SF", "FSTL_1_D_105 SF", "FSTL_1_D_105 750nM SF", 
                     "FSTL_1_E_097 100nM SF", "FSTL_1_E_097 SF", "FSTL_1_E_097 750nM SF", "FSTL_1_F_109 100nM SF", "FSTL_1_F_109 SF", "FSTL_1_F_109 750nM SF", 
                     "FSTL_1_H_201 100nM SF", "FSTL_1_H_201 SF", "FSTL_1_H_201 750nM SF", "FSTL_1_I_204 100nM SF", "FSTL_1_I_204 SF", "FSTL_1_I_204 750nM SF", 
                     "FSTL_1_K_219 100nM SF", "FSTL_1_K_219 SF", "FSTL_1_K_219 750nM SF", "FSTL_1_*_259 100nM SF", "FSTL_1_*_259 SF", "FSTL_1_*_259 750nM SF", 
                     "FSTL_1_L_412 100nM SF", "FSTL_1_L_412 SF", "FSTL_1_L_412 750nM SF", "FSTL_1_M_398 100nM SF", "FSTL_1_M_398 SF", "FSTL_1_M_398 750nM SF", 
                     "FSTL_1_P_820 100nM SF", "FSTL_1_P_820 SF", "FSTL_1_P_820 750nM SF", "FSTL_1_T_770 100nM SF", "FSTL_1_T_770 SF", "FSTL_1_T_770 750nM SF", 
                     "FSTL_1_F_109 100nM LC", "FSTL_1_F_109 LC", "FSTL_1_F_109 750nM LC", "FSTL_1_H_201 100nM LC", "FSTL_1_H_201 LC", "FSTL_1_H_201 750nM LC", 
                     "FSTL_1_I_204 100nM LC", "FSTL_1_I_204 LC", "FSTL_1_I_204 750nM LC", "FSTL_1_K_219 100nM LC", "FSTL_1_K_219 LC", "FSTL_1_K_219 750nM LC", 
                     "FSTL_1_M_398 100nM LC", "FSTL_1_M_398 LC", "FSTL_1_M_398 750nM LC", "FSTL_1_O_417 100nM LC", "FSTL_1_O_417 LC", "FSTL_1_O_417 750nM LC", 
                     "FSTL_1_O_417 100nM SF", "FSTL_1_O_417 SF", "FSTL_1_O_417 750nM SF", "FSTL_1_P_820 100nM LC", "FSTL_1_P_820 LC", "FSTL_1_P_820 750nM LC"]
        


        if True:
            saveTD0 = {}
            saveTD0fluor = {}
            eff_perTarget = {}


            pTab = [["run", header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar], header[0][rar_excl], "start", "end", "TD0", "TD0 fluor", "TD0 end", "4cycInc", "BG", "Baseline Err", "Dilution Eff", "indiv PCR Eff"]]
            for row in range(0, spFl[0]):
                if res[row][rar_sample_type] in ["ntc", "nac", "ntp", "nrt", "opt"]:
                    continue

                dilStd = ""
                if res[row][rar_tar] in effDiutions:
                    dilStd = str(effDiutions[res[row][rar_tar]])
                pTab.append([self["id"] , res[row][rar_id], res[row][rar_well], res[row][rar_sample], res[row][rar_tar], res[row][rar_excl], rawMinBySD[row], maxSD[row], rawTD0[row], 
                             td0_fluor[row], int(rawTD0[row]), check_four_cyc_inc_TD0[row], backgroundNew[row], baselineError[row], dilStd, indiv_PCR_Eff[row]])
                if baselineError[row] == 0:
                    if res[row][rar_tar] not in eff_perTarget:
                        eff_perTarget[res[row][rar_tar]] = {}
                        eff_perTarget[res[row][rar_tar]]["raw"] = []
                        eff_perTarget[res[row][rar_tar]]["r2"] = []
                        eff_perTarget[res[row][rar_tar]]["eff"] = []
                        eff_perTarget[res[row][rar_tar]]["mean_eff"] = []
                        eff_perTarget[res[row][rar_tar]]["wol_eff"] = []
                        eff_perTarget[res[row][rar_tar]]["excl_indiv"] = []
                        eff_perTarget[res[row][rar_tar]]["excl_wol"] = []
                        eff_perTarget[res[row][rar_tar]]["excl_sample"] = []
                        eff_perTarget[res[row][rar_tar]]["excl_eff_indiv"] = []
                        eff_perTarget[res[row][rar_tar]]["excl_eff_wol"] = []
                        eff_perTarget[res[row][rar_tar]]["diff_mean_eff"] = []
                        eff_perTarget[res[row][rar_tar]]["sum_diff_eff"] = []
                        eff_perTarget[res[row][rar_tar]]["well"] = []
                        eff_perTarget[res[row][rar_tar]]["sample"] = []

                    eff_perTarget[res[row][rar_tar]]["well"].append(res[row][rar_well])
                    eff_perTarget[res[row][rar_tar]]["sample"].append(res[row][rar_sample])


                    eff_perTarget[res[row][rar_tar]]["raw"].append(indiv_PCR_Eff[row])

                    stop = int(rawTD0[row])
                    # At minimum 4 cycles are required
                    if stop < 4:
                        baselineError[row] = 1
                        continue
                    # At maximum 9 cycles should be inluded
                    start = max(rawMinBySD[row] + 1, stop - 9 + 1)
                    # Now check that 4 cycles are left
                    if stop - start + 1 < 4:
                        baselineError[row] = 2
                        continue
                    eff_arr = []

                    slope, intercept, r_value, p_value, std_err = scp.stats.linregress(range(start - 1, stop), baselineCorrectedData[row][start - 1:stop])
                    eff_perTarget[res[row][rar_tar]]["r2"].append(r_value)

                    for count in range(0, stop - start - 2):
                        selSection = baselineCorrectedData[row][stop - count - 4 : stop -count]
                        eff = np.power(10, _lrp_testSlopes3(selSection))
                        eff_arr.append(eff)
                    meaneff = np.mean(eff_arr)
                    diffmean = eff_arr - meaneff
                    sumDiff = np.sum(np.abs(diffmean)) / len(diffmean)

                    eff_perTarget[res[row][rar_tar]]["eff"].append(eff_arr)
                    eff_perTarget[res[row][rar_tar]]["mean_eff"].append(meaneff)
                    eff_perTarget[res[row][rar_tar]]["wol_eff"].append(WoL_PCR_Eff[row])
                    eff_perTarget[res[row][rar_tar]]["excl_indiv"].append(vecExclIndivEff[row])
                    eff_perTarget[res[row][rar_tar]]["excl_wol"].append(vecExclEff[row])
                    eff_perTarget[res[row][rar_tar]]["excl_sample"].append(vecSkipSample[row])
                    eff_perTarget[res[row][rar_tar]]["excl_eff_indiv"].append(mean_indiv_PCR_Eff[row])
                    eff_perTarget[res[row][rar_tar]]["excl_eff_wol"].append(mean_PCR_Eff[row])
                    eff_perTarget[res[row][rar_tar]]["diff_mean_eff"].append(diffmean)
                    eff_perTarget[res[row][rar_tar]]["sum_diff_eff"].append(sumDiff)

                    if res[row][rar_tar] not in saveTD0:
                        saveTD0[res[row][rar_tar]] = {}
                        saveTD0[res[row][rar_tar]]["raw"] = []
                        saveTD0fluor[res[row][rar_tar]] = {}
                        saveTD0fluor[res[row][rar_tar]]["raw"] = []
                    saveTD0[res[row][rar_tar]]["raw"].append(rawTD0[row])
                    saveTD0fluor[res[row][rar_tar]]["raw"].append(td0_fluor[row])

          #      for bla in dddddd[row]:
           #         pTab[row + 1].append(str(bla))  vecExclIndivEff
                
            
            ww = open("/home/untergasser/code/rdml/tools/server/rdmlpython/test/temp_td0.csv", "w")
            for row in range(0, len(pTab)):
                for col in range(0, len(pTab[row])):
                    ww.write(str(pTab[row][col]) + "\t")
                ww.write("\n")
            ww.close()

            ww = open("/home/untergasser/code/rdml/tools/server/rdmlpython/test/temp_efficency.csv", "w")
            ww.write("Target\t" + "Sample\t" + "Well\t" + "Dilution Efficiency\t" + "WoL Efficiency\t" + "Mean Curve Efficiency\t" + "Out Indiv Eff\t" +
                     "Out WoL Eff\t" + "Outlier Individ\t" + "Outlier WoL\t" + "Skip Sample\t" + "R^2\t" + "Individ Eff\t" + "\t" + "Mean Window Eff\t" + "\t")
            for col in range(0, 7):
                ww.write(str(col) + " - " + str(col + 4))
                ww.write("\t")
            ww.write("\n")
            for tar in eff_orddr:
                if tar not in eff_perTarget:
                    continue
                for pos in range(0, len(eff_perTarget[tar]["sample"])):
                    ww.write(tar + "\t")
                    ww.write(str(eff_perTarget[tar]["sample"][pos]) + "\t")
                    ww.write(str(eff_perTarget[tar]["well"][pos]) + "\t")
                    if tar in effDiutions:
                        ww.write(str(effDiutions[tar]))
                    ww.write("\t")
                    ww.write("{:.4f}".format(np.mean(eff_perTarget[tar]["wol_eff"])) + "\t")
                    ww.write("{:.4f}".format(np.mean(eff_perTarget[tar]["raw"])) + "\t")
                    ww.write("{:.4f}".format(eff_perTarget[tar]["excl_eff_indiv"][pos]) + "\t")
                    ww.write("{:.4f}".format(eff_perTarget[tar]["excl_eff_wol"][pos]) + "\t")
                    ww.write(str(eff_perTarget[tar]["excl_indiv"][pos]) + "\t")
                    ww.write(str(eff_perTarget[tar]["excl_wol"][pos]) + "\t")
                    ww.write(str(eff_perTarget[tar]["excl_sample"][pos]) + "\t")
                    ww.write("{:.6f}".format(eff_perTarget[tar]["r2"][pos]) + "\t")
                    ww.write("{:.4f}".format(eff_perTarget[tar]["raw"][pos]) + "\t")
                    ww.write("\t")
                    ww.write("{:.4f}".format(np.mean(eff_perTarget[tar]["mean_eff"][pos])) + "\t")
                    ww.write("\t")
                    for col in range(0, 7):
                        if col < len(eff_perTarget[tar]["eff"][pos]):
                            ww.write("{:.4f}".format(eff_perTarget[tar]["eff"][pos][col]))
                        ww.write("\t")
                    ww.write("\n")
                ww.write("\n")
                ww.write("\n")
            ww.close()

            for tar in eff_perTarget:
                eff_perTarget[tar]["mean"] = np.mean(eff_perTarget[tar]["raw"])
                eff_perTarget[tar]["std"] = np.std(eff_perTarget[tar]["raw"], ddof=1)
            for tar in saveTD0:
                saveTD0[tar]["mean"] = np.mean(saveTD0[tar]["raw"])
                saveTD0[tar]["std"] = np.std(saveTD0[tar]["raw"], ddof=1)
                saveTD0fluor[tar]["mean"] = np.mean(saveTD0fluor[tar]["raw"])
                saveTD0fluor[tar]["std"] = np.std(saveTD0fluor[tar]["raw"], ddof=1)

            if 0:
                for tar in saveTD0:
                    if tar in saveTD0:
                        print(tar + " " + "{:.4f}".format(saveTD0[tar]["mean"]) + " (" + "{:.4f}".format(saveTD0[tar]["std"]) + ") " + "{:.4f}".format(saveTD0fluor[tar]["mean"]) + " (" + "{:.4f}".format(saveTD0fluor[tar]["std"]) + ")")


            for tar in eff_orddr:
                dil = -1.0
                if tar in effDiutions:
                    dil = effDiutions[tar]
                if tar in eff_perTarget:
                    if 0:
                        print(tar + " " + "{:.4f}".format(dil) + " " + "{:.4f}".format(dil - eff_perTarget[tar]["mean"]) + "   " + "{:.4f}".format(eff_perTarget[tar]["mean"]) + " (" + "{:.4f}".format(eff_perTarget[tar]["std"]) + ")")


         #   print("--------------------------")







        #########################
        # write out the results #
        #########################
        for rRow in range(0, len(res)):
            res[rRow][rar_baseline] = vecBackground[rRow] - negShiftBaseline[rRow]
            res[rRow][rar_plat_val] = meanPlateau[rRow]
            if res[rRow][rar_baseline] != 0.0:
                res[rRow][rar_plat_base] = meanPlateau[rRow] / res[rRow][rar_baseline]
            else:
                res[rRow][rar_plat_base] = 0.0
            res[rRow][rar_lower_limit] = lowWin[vecTarget[rRow]]
            res[rRow][rar_upper_limit] = upWin[vecTarget[rRow]]

            res[rRow][rar_n_log] = stopCyc[rRow] - startCycFix[rRow] + 1
            res[rRow][rar_stop_log] = stopCyc[rRow]
            res[rRow][rar_n_included] = nInclu[rRow]

            res[rRow][rar_indiv_PCR_eff] = indiv_PCR_Eff[rRow]
            res[rRow][rar_R2] = correl[rRow] * correl[rRow]
            res[rRow][rar_PCR_eff] = mean_PCR_Eff[rRow]
            res[rRow][rar_PCR_eff_err] = mean_PCR_Eff_Err[rRow]

            res[rRow][rar_threshold] = indiv_thres[rRow]
            res[rRow][rar_Cq] = rawTD0[rRow]
            res[rRow][rar_indiv_Ncopy] = -1.0
            res[rRow][rar_Ncopy] = -1.0
            ss_fix = 1.0
            if res[rRow][rar_sample_nucleotide] == "ss":
                ss_fix = 2.0
            res[rRow][rar_nAmpli] = ss_fix * target_limit[res[rRow][rar_tar]] * reactVol[rRow]
            if rawTD0[rRow] > 0.0:
                if 2.3 > WoL_PCR_Eff[rRow] > 1.0:
                    res[rRow][rar_indiv_Ncopy] = res[rRow][rar_nAmpli] / np.power(WoL_PCR_Eff[rRow], rawTD0[rRow])
                if 2.3 > mean_PCR_Eff[rRow] > 1.0:
                    res[rRow][rar_Ncopy] = res[rRow][rar_nAmpli] / np.power(mean_PCR_Eff[rRow], rawTD0[rRow])

         #   if td0_Cq[rRow] > 0.0:
          #      print(str(res[rRow][1]) + " - " +  str(res[rRow][rar_Ncopy])  + " - " +  str(td0_Cq[rRow]) + " - " +  str(mean_PCR_Eff[rRow]))
 
            res[rRow][rar_amplification] = not vecNoAmplification[rRow]
            res[rRow][rar_baseline_error] = vecBaselineError[rRow]
            res[rRow][rar_instable_baseline] = vecInstableBaseline[rRow]
            res[rRow][rar_plateau] = not vecNoPlateau[rRow]
            res[rRow][rar_noisy_sample] = vecNoisySample[rRow]
            res[rRow][rar_excl_eff] = vecExclEff[rRow]
            res[rRow][rar_shortLogLinPhase] = vecShortLogLin[rRow]
            res[rRow][rar_CqIsShifting] = vecCtIsShifting[rRow]
            res[rRow][rar_tooLowCqEff] = vecTooLowCqEff[rRow]
            res[rRow][rar_tooLowCqN0] = vecTooLowCqN0[rRow]
            res[rRow][rar_isUsedInWoL] = vecIsUsedInWoL[rRow]

        # Delete
            res[rRow][del_threshold] = threshold[0]
            res[rRow][del_N0] = nNull[rRow]

        mmN0_sum = 0.0
        mmNcopy_sum = 0.0
        mmN0_num = 0
        mmNcopy_num = 0
        geomeanNo = 1.0
        geomeanNcopy = 1.0
        for rRow in range(0, len(res)):
            if res[rRow][rar_isUsedInWoL]:
                if res[rRow][rar_Ncopy] > 0.0 and res[rRow][del_N0] > 0.0:
                    mmN0_sum += math.log(res[rRow][del_N0])
                    mmNcopy_sum += math.log(res[rRow][rar_Ncopy])
                    mmN0_num += 1
                    mmNcopy_num += 1
        if mmNcopy_num > 0:
            geomeanNo = math.exp(mmN0_sum / mmN0_num)
            geomeanNcopy = math.exp(mmNcopy_sum / mmNcopy_num)

        for rRow in range(0, len(res)):
            if nNull[rRow] > 0.0:
                res[rRow][del_ncopy] = nNull[rRow] * geomeanNcopy / geomeanNo
            else:
                res[rRow][del_ncopy] = -1.0

            if res[rRow][rar_Ncopy] == res[rRow][del_ncopy]:
                    res[rRow][del_fact] = 1.0
            elif res[rRow][rar_Ncopy] <= 0.0:
                res[rRow][del_fact] = -0.1
            elif res[rRow][del_ncopy] <= 0.0:
                res[rRow][del_fact] = -0.2
            else:
                if res[rRow][del_ncopy] < res[rRow][rar_Ncopy]:
                    res[rRow][del_fact] = res[rRow][del_ncopy] / res[rRow][rar_Ncopy]
                else:
                    res[rRow][del_fact] = res[rRow][rar_Ncopy] / res[rRow][del_ncopy]


        # Jans Columns
        header[0].append("0: mean efficiency")
        header[0].append("1: baseline value")
        header[0].append("2: individual efficiency")
        header[0].append("3: threshold value (as set by RDML 1")
        header[0].append("4: Cq value as crossing point between fluorescence curve and threshold (3)")
        header[0].append("5: Cq value as crossing point between the line fitted to the curve with individual efficiency and threshold (3)")
        header[0].append("6: Cq value as crossing point between the line fitted to the curve with the average efficiency and threshold (3)")
        header[0].append("7: TD0 calculated as the zero-crossing in the third derivative that is closest to the maximum in the second derivative")
        header[0].append("8: Fluorescence at TD0, calculated by interpolation on the observed fluorescence curve")
        header[0].append("9: Threshold calculated as average fluorescence at TD0 (column 8)")
        header[0].append("10: TD0 value as crossing point between fluorescence curve and threshold 9")
        header[0].append("11: TD0 value as crossing point between the line fitted with individual efficiency to the curve and threshold 9")
        header[0].append("12: TD0 value as crossing point between the line fitted with the average efficiency to the curve and  threshold 9")
        header[0].append("13: average fluorescence values at plateau level (5 cycles, excluding the last cycle)")
        header[0].append("14: standard deviation of plateau fluorescence")
        header[0].append("15: N0 based on individual efficiency (column 2), threshold (column 3) and Cq from column 4")
        header[0].append("16: N0 based on individual efficiency (column 2), threshold (column 3) and Cq from column 5")
        header[0].append("17: N0 based on individual efficiency (column 2), threshold (column 3) and Cq from column 6")
        header[0].append("18: N0 based on average efficiency (0), threshold and Cq from column 4")
        header[0].append("19: N0 based on average efficiency (0), threshold and Cq from column 5")
        header[0].append("20: N0 based on average efficiency (0), threshold and Cq from column 6")
        header[0].append("21: Ncopy based on TD0 from column 7 and individual efficiency (column 2)")
        header[0].append("22: Ncopy based on TD0 from column 10 and individual efficiency (column 2)")
        header[0].append("23: Ncopy based on TD0 from column 11 and individual efficiency (column 2)")
        header[0].append("24: Ncopy based on TD0 from column 12 and individual efficiency (column 2)")
        header[0].append("25: Ncopy based on TD0 from column 7 and average efficiency (column 0)")
        header[0].append("26: Ncopy based on TD0 from column 10 and average efficiency (column 0)")
        header[0].append("27: Ncopy based on TD0 from column 11 and average efficiency (column 0)")
        header[0].append("28: Ncopy based on TD0 from column 12 and average efficiency (column 0)")


        for rRow in range(0, len(res)):
            res[rRow].append(mean_PCR_Eff[rRow])  # 0
            res[rRow].append(res[rRow][rar_baseline])  # 1
            res[rRow].append(res[rRow][rar_indiv_PCR_eff])  # 2
            res[rRow].append(threshold[0])  # 3
            xLast = -1
            xOver = 0
            for cyci in range(1, len(baselineCorrectedData[rRow])):
                xLast += 1
                xOver += 1
                if baselineCorrectedData[rRow, cyci] > threshold[0]:
                    break
            xxInc = np.log10(baselineCorrectedData[rRow, xOver]) - np.log10(baselineCorrectedData[rRow, xLast])
            xxBase = np.log10(threshold[0]) - np.log10(baselineCorrectedData[rRow, xLast])
            if xxInc > 0.0:
                xxCq = 1.0 + xLast + xxBase / xxInc
            else:
                xxCq = 1.0 + xLast
            res[rRow].append(xxCq)  # 4
            res[rRow].append(indivCq[rRow])  # 5
            res[rRow].append(meanCq[rRow])  # 6
            res[rRow].append(rawTD0[rRow])  # 7
            res[rRow].append(td0_fluor[rRow])  # 8
            res[rRow].append(td0_mean_thres[rRow])  # 9
            xLast = -1
            xOver = 0
            for cyci in range(1, len(baselineCorrectedData[rRow])):
                xLast += 1
                xOver += 1
                if baselineCorrectedData[rRow, cyci] > td0_mean_thres[rRow]:
                    break
            xxInc = np.log10(baselineCorrectedData[rRow, xOver]) - np.log10(baselineCorrectedData[rRow, xLast])
            if td0_mean_thres[rRow] > 0.0 and baselineCorrectedData[rRow, xLast] > 0.0:
                xxBase = np.log10(td0_mean_thres[rRow]) - np.log10(baselineCorrectedData[rRow, xLast])
            else:
                xxBase = -1.0
            if xxInc > 0.0:
                xxTD0 = 1.0 + xLast + xxBase / xxInc
            else:
                xxTD0 = 1.0 + xLast
            res[rRow].append(xxTD0)  # 10
            res[rRow].append(td0_indiCq[rRow])  # 11
            if td0_mean_thres[rRow] > 0.0 and mean_PCR_Eff[rRow] > 0.0:
                td0_meanCq[rRow] = indMeanX[rRow] + (np.log10(td0_mean_thres[rRow]) - indMeanY[rRow]) / np.log10(mean_PCR_Eff[rRow])
            else:
                td0_meanCq[rRow] = -1.0
            res[rRow].append(td0_meanCq[rRow])  # 12
            xxFluorLen = len(baselineCorrectedData[rRow])
            xxPlatLast = baselineCorrectedData[rRow, xxFluorLen - 6 : xxFluorLen - 1]
            res[rRow].append(np.mean(xxPlatLast))  # 13
            res[rRow].append(np.std(xxPlatLast))  # 14
            if res[rRow][rar_indiv_PCR_eff] > 0.0 and xxCq > 0.0:
                xxnNull = threshold[0] / np.power(res[rRow][rar_indiv_PCR_eff], xxCq)
                res[rRow].append(xxnNull)  # 15
            else:
                res[rRow].append(-1.0)
            if res[rRow][rar_indiv_PCR_eff] > 0.0 and indivCq[rRow] > 0.0:
                xxnNull = threshold[0] / np.power(res[rRow][rar_indiv_PCR_eff], indivCq[rRow])
                res[rRow].append(xxnNull)  # 16
            else:
                res[rRow].append(-1.0)
            if res[rRow][rar_indiv_PCR_eff] > 0.0 and meanCq[rRow] > 0.0:
                xxnNull = threshold[0] / np.power(res[rRow][rar_indiv_PCR_eff], meanCq[rRow])
                res[rRow].append(xxnNull)  # 17
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and 100 > xxCq > 0.0:
                xxnNull = threshold[0] / np.power(mean_PCR_Eff[rRow], xxCq)
                res[rRow].append(xxnNull)  # 18
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and 100 > indivCq[rRow] > 0.0:
                xxnNull = threshold[0] / np.power(mean_PCR_Eff[rRow], indivCq[rRow])
                res[rRow].append(xxnNull)  # 19
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and 100 > meanCq[rRow] > 0.0:
                xxnNull = threshold[0] / np.power(mean_PCR_Eff[rRow], meanCq[rRow])
                res[rRow].append(xxnNull)  # 20
            else:
                res[rRow].append(-1.0)

            if WoL_PCR_Eff[rRow] > 0.0 and rawTD0[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(WoL_PCR_Eff[rRow], rawTD0[rRow])
                res[rRow].append(xxCopy)  # 21
            else:
                res[rRow].append(-1.0)
            if WoL_PCR_Eff[rRow] > 0.0 and 100 > xxTD0 > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(WoL_PCR_Eff[rRow], xxTD0)
                res[rRow].append(xxCopy)  # 22
            else:
                res[rRow].append(-1.0)
            if WoL_PCR_Eff[rRow] > 0.0 and td0_indiCq[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(WoL_PCR_Eff[rRow], td0_indiCq[rRow])
                res[rRow].append(xxCopy)  # 23
            else:
                res[rRow].append(-1.0)
            if WoL_PCR_Eff[rRow] > 0.0 and td0_meanCq[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(WoL_PCR_Eff[rRow], td0_meanCq[rRow])
                res[rRow].append(xxCopy)  # 24
            else:
                res[rRow].append(-1.0)

            if mean_PCR_Eff[rRow] > 0.0 and rawTD0[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(mean_PCR_Eff[rRow], rawTD0[rRow])
                res[rRow].append(xxCopy)  # 25
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and 100 > xxTD0 > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(mean_PCR_Eff[rRow], xxTD0)
                res[rRow].append(xxCopy)  # 26
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and td0_indiCq[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(mean_PCR_Eff[rRow], td0_indiCq[rRow])
                res[rRow].append(xxCopy)  # 27
            else:
                res[rRow].append(-1.0)
            if mean_PCR_Eff[rRow] > 0.0 and td0_meanCq[rRow] > 0.0:
                xxCopy = res[rRow][rar_nAmpli] / np.power(mean_PCR_Eff[rRow], td0_meanCq[rRow])
                res[rRow].append(xxCopy)  # 28
            else:
                res[rRow].append(-1.0)


        ###################################
        # calculate excl and note strings #
        ###################################
        for rRow in range(0, len(res)):
            exclVal = _cleanErrorString(res[rRow][rar_excl], "amp")
            noteVal = _cleanErrorString(res[rRow][rar_note], "amp")

            ncopy = res[rRow][rar_Ncopy]
            diffMeanEff = False #  res[rRow][rar_effOutlier_Skip_Plat_Out]

            if res[rRow][rar_sample_type] in ["ntc", "nac", "ntp", "nrt"]:
                if ncopy > 0.0:
                    exclVal += "amplification in negative control;"
                    if res[rRow][rar_baseline_error]:
                        exclVal += "baseline error;"
                else:
                    if res[rRow][rar_amplification]:
                        exclVal += "amplification in negative control;"
                    if res[rRow][rar_baseline_error]:
                        exclVal += "baseline error;"
                if res[rRow][rar_plateau]:
                    noteVal += "plateau in negative control;"

            if res[rRow][rar_sample_type] in ["std", "pos"]:
                if not (ncopy > 0.0):
                    exclVal += "no amplification in positive control;"
                else:
                    if not res[rRow][rar_amplification]:
                        exclVal += "no amplification in positive control;"
                if res[rRow][rar_baseline_error]:
                    exclVal += "baseline error in positive control;"
                if res[rRow][rar_instable_baseline]:
                    exclVal += "instable baseline in positive control;"
                if not res[rRow][rar_plateau]:
                    noteVal += "no plateau in positive control;"
                if res[rRow][rar_noisy_sample]:
                    noteVal += "noisy sample in positive control;"

                if ncopy > 100000000.0:
                    noteVal += "Ncopy > 100 Mio;"
                if ncopy < 10.0:
                    noteVal += "Ncopy < 10;"
                if res[rRow][rar_n_log] < 5:
                    noteVal += "only " + str(res[rRow][rar_n_log]) + " values in log phase;"
                if res[rRow][rar_indiv_PCR_eff] < 1.7:
                    noteVal += "indiv PCR eff is " + "{:.3f}".format(res[rRow][rar_indiv_PCR_eff]) + " < 1.7;"
                if diffMeanEff:
                    if not np.isfinite(res[rRow][rar_indiv_PCR_eff]):
                        noteVal += "no indiv PCR eff can be calculated;"
                    else:
                        if excludeEfficiency == "outlier":
                            noteVal += "PCR efficiency outlier;"
                        diffFromMean = res[rRow][rar_indiv_PCR_eff] - res[rRow][rar_PCR_eff]
                        if diffFromMean > 0.0:
                            noteVal += "indiv PCR eff is higher than mean PCR eff by "
                            noteVal += "{:.3f}".format(diffFromMean) + ";"
                        else:
                            noteVal += "indiv PCR eff is lower than mean PCR eff by "
                            noteVal += "{:.3f}".format(-1 * diffFromMean) + ";"

            if res[rRow][rar_sample_type] in ["unkn"]:
                if not res[rRow][rar_amplification]:
                    noteVal += "no amplification;"
                if res[rRow][rar_baseline_error]:
                    noteVal += "baseline error;"
                if res[rRow][rar_instable_baseline]:
                    noteVal += "instable baseline;"
                if not res[rRow][rar_plateau]:
                    noteVal += "no plateau;"
                if res[rRow][rar_noisy_sample]:
                    noteVal += "noisy sample;"

                if ncopy > 100000000.0:
                    noteVal += "Ncopy > 100 Mio;"
                if 0.0 < ncopy < 10.0:
                    noteVal += "Ncopy < 10;"
                if res[rRow][rar_n_log] < 5:
                    noteVal += "only " + str(res[rRow][rar_n_log]) + " values in log phase;"
                if res[rRow][rar_indiv_PCR_eff] < 1.7:
                    noteVal += "indiv PCR eff is " + "{:.3f}".format(res[rRow][rar_indiv_PCR_eff]) + " < 1.7;"
                if diffMeanEff:
                    if not np.isfinite(res[rRow][rar_indiv_PCR_eff]):
                        noteVal += "no indiv PCR eff can be calculated;"
                    else:
                        if excludeEfficiency == "outlier":
                            noteVal += "PCR efficiency outlier;"
                        diffFromMean = res[rRow][rar_indiv_PCR_eff] - meanEffVal
                        if diffFromMean > 0.0:
                            noteVal += "indiv PCR eff is higher than mean PCR eff by "
                            noteVal += "{:.3f}".format(diffFromMean) + ";"
                        else:
                            noteVal += "indiv PCR eff is lower than mean PCR eff by "
                            noteVal += "{:.3f}".format(-1 * diffFromMean) + ";"

            # Write back
            exclVal = re.sub(r'^;|;$', '', exclVal)
            noteVal = re.sub(r'^;|;$', '', noteVal)
            res[rRow][rar_excl] = exclVal
            res[rRow][rar_note] = noteVal

        ##############################
        # write out the rdml results #
        ##############################
        # updateRDML = True
        menJanTres = np.mean(indiv_thres)
        if updateRDML is True:
            self["backgroundDeterminationMethod"] = 'LinRegPCR, constant'
            self["cqDetectionMethod"] = 'automated threshold and baseline settings'
            dataXMLelements = _getXMLDataType()
            collectedTargetEff = {}
            collectedTargetErr = {}
            for rRow in range(0, len(res)):
                if rdmlElemData[rRow] is not None:
                    cqVal = rawTD0[rRow]
                    meanEffVal = mean_PCR_Eff[rRow]
                    stdEffVal = mean_PCR_Eff_Err[rRow]
                    if np.isnan(cqVal) or cqVal < 0.0 or cqVal > 1000.0:
                        goodVal = "-1.0"
                    else:
                        goodVal = "{:.3f}".format(cqVal)

                    _change_subelement(rdmlElemData[rRow], "cq", dataXMLelements, goodVal, True, "string")
                    _change_subelement(rdmlElemData[rRow], "excl", dataXMLelements, res[rRow][rar_excl], True, "string")
                    if dataVersion in ["1.3", "1.4"]:
                        _change_subelement(rdmlElemData[rRow], "N0", dataXMLelements, "", True, "string")
                        if dataVersion in ["1.4"]:
                            goodVal = res[rRow][rar_Ncopy]
                            if np.isfinite(goodVal):
                                if goodVal > 0.0:
                                    goodVal = "{:.4f}".format(goodVal)
                                    _change_subelement(rdmlElemData[rRow], "Ncopy", dataXMLelements, goodVal, True, "string")
                                else:
                                    _change_subelement(rdmlElemData[rRow], "Ncopy", dataXMLelements, "-1.0", True, "string")
                        _change_subelement(rdmlElemData[rRow], "ampEffMet", dataXMLelements, "LinRegPCR", True, "string")
                        if not np.isfinite(meanEffVal):
                            goodVal = "-1.0"
                        else:
                            goodVal = "{:.3f}".format(meanEffVal)
                        _change_subelement(rdmlElemData[rRow], "ampEff", dataXMLelements, goodVal, True, "string")
                        if updateTargetEfficiency:
                            collectedTargetEff[res[rRow][rar_tar]] = goodVal
                        if not np.isfinite(stdEffVal):
                            goodVal = "-1.0"
                        else:
                            goodVal = "{:.3f}".format(stdEffVal)
                        _change_subelement(rdmlElemData[rRow], "ampEffSE", dataXMLelements, goodVal, True, "string")
                        if updateTargetEfficiency:
                             collectedTargetErr[res[rRow][rar_tar]] = goodVal
                        _change_subelement(rdmlElemData[rRow], "note", dataXMLelements, res[rRow][rar_note], True, "string")
                    goodVal = "{:.3f}".format(vecBackground[rRow] - negShiftBaseline[rRow])
                    _change_subelement(rdmlElemData[rRow], "bgFluor", dataXMLelements, goodVal, True, "string")
                    goodVal = "{:.3f}".format(threshold[0])
                    _change_subelement(rdmlElemData[rRow], "quantFluor", dataXMLelements, goodVal, True, "string")
            if updateTargetEfficiency:
                tarXMLKeys = ["description", "documentation", "xRef", "type", "amplificationEfficiencyMethod",
                              "amplificationEfficiency", "amplificationEfficiencySE", "meltingTemperature",
                              "detectionLimit", "dyeId", "sequences", "commercialAssay"]
                for curTar in collectedTargetEff:
                    eleTars = _get_all_children(rootPar, "target")
                    for eleTar in eleTars:
                        foundId = eleTar.attrib['id']
                        if foundId != curTar:
                            continue
                        eleEff = _get_or_create_subelement(eleTar, "amplificationEfficiency", tarXMLKeys)
                        eleEff.text = collectedTargetEff[curTar]
                        if dataVersion != "1.1":
                            if curTar in collectedTargetErr:
                                eleEffSE = _get_or_create_subelement(eleTar, "amplificationEfficiencySE", tarXMLKeys)
                                eleEffSE.text = collectedTargetErr[curTar]

        if timeRun:
            stop_time = datetime.datetime.now() - start_time
            print("Done All: " + str(stop_time) + "sec")

        if saveResultsCSV:
            retCSV = ""
            for rCol in range(0, len(header[0])):
                retCSV += header[0][rCol] + "\t"
            retCSV = re.sub(r"\t$", "\n", retCSV)

            for rRow in range(0, len(res)):
                for rCol in range(0, len(res[rRow])):
                    if rCol in [rar_amplification, rar_baseline_error, rar_instable_baseline, rar_plateau, rar_noisy_sample,
                                rar_excl_eff, rar_shortLogLinPhase, rar_CqIsShifting,
                                rar_tooLowCqEff, rar_tooLowCqN0, rar_isUsedInWoL]:
                        if res[rRow][rCol]:
                            retCSV += "Yes\t"
                        else:
                            retCSV += "No\t"
                    elif rCol in [rar_baseline, rar_lower_limit, rar_upper_limit,
                                  rar_indiv_PCR_eff, rar_R2, rar_PCR_eff, rar_PCR_eff_err]:
                        retCSV += "{0:0.6f}".format(float(res[rRow][rCol])) + "\t"
                    elif rCol in [rar_Cq, rar_indiv_Ncopy, rar_Ncopy]:
                        retCSV += "{0:0.4f}".format(float(res[rRow][rCol])) + "\t"
                      # retCSV += "{0:0.6e}".format(float(res[rRow][rCol])) + "\t"
                    else:
                        retCSV += str(res[rRow][rCol]) + "\t"
                retCSV = re.sub(r"\t$", "\n", retCSV)
            finalData["resultsCSV"] = retCSV

        if saveResultsList:
            finalData["resultsList"] = header + res

        if err != "":
            finalData["error"] = err

        return finalData

    def webAppMeltCurveAnalysis(self, normMethod="exponential", fluorSource="normalised",
                                truePeakWidth=1.0, artifactPeakWidth=1.0,
                                expoLowTemp=65.0, expoHighTemp=92.0,
                                bilinLowStartTemp=65.0, bilinLowStopTemp=67.0,
                                bilinHighStartTemp=93.0, bilinHighStopTemp=94.0,
                                peakLowTemp=60.0, peakHighTemp=98.0,
                                peakMaxWidth=5.0, peakCutoff=5.0,
                                updateRDML=False):
        """Performs LinRegPCR on the run. Modifies the cq values and returns a json with additional data.

        Args:
            self: The class self parameter.
            normMethod: The normalization method "exponential", "bilinear" or "combined"
            fluorSource: choose "normalised" or "smooth" for fluorescence correction factor calculation
            truePeakWidth: the maximum allowed width in Celsius of the expected (true) peak
            artifactPeakWidth: the maximum allowed width in Celsius of all artifact peaks
            expoLowTemp: the low temperature for the exponential normalisation
            expoHighTemp: the high temperature for the exponential normalisation
            bilinLowStartTemp: the low start temperature for the bilinear normalisation
            bilinLowStopTemp: the low stop temperature for the bilinear normalisation
            bilinHighStartTemp: the high start temperature for the bilinear normalisation
            bilinHighStopTemp: the high stop temperature for the bilinear normalisation
            peakLowTemp: peaks below this temperature are ignored
            peakHighTemp: peaks above this temperature are ignored
            peakMaxWidth: peaks broader than this temperature width are ignored
            peakCutoff: the percentage below melting peaks are ignored in calculations
            updateRDML: If true, update the RDML data with the calculated values.

        Returns:
            A dictionary with the resulting data, presence and format depending on input.
            rawData: A 2d array with the raw fluorescence values
            baselineCorrectedData: A 2d array with the baseline corrected raw fluorescence values
            resultsList: A 2d array object.
            resultsCSV: A csv string.
        """

        allData = self.getreactjson()
        res = self.meltCurveAnalysis(normMethod=normMethod, fluorSource=fluorSource,
                                     truePeakWidth=truePeakWidth, artifactPeakWidth=artifactPeakWidth,
                                     expoLowTemp=expoLowTemp, expoHighTemp=expoHighTemp,
                                     bilinLowStartTemp=bilinLowStartTemp, bilinLowStopTemp=bilinLowStopTemp,
                                     bilinHighStartTemp=bilinHighStartTemp, bilinHighStopTemp=bilinHighStopTemp,
                                     peakLowTemp=peakLowTemp, peakHighTemp=peakHighTemp,
                                     peakMaxWidth=peakMaxWidth, peakCutoff=peakCutoff, updateRDML=updateRDML,
                                     saveRaw=False, saveDerivative=True,
                                     saveResultsList=True, saveResultsCSV=False, verbose=False)
        if "derivative" in res:
            if "smoothed" in res["derivative"]:
                bas_temp_min = 120.0
                bas_temp_max = 0.0
                bas_fluor_min = 99999999
                bas_fluor_max = 0.0
                for row in range(1, len(res["derivative"]["smoothed"])):
                    bass_json = []
                    for col in range(6, len(res["derivative"]["smoothed"][row])):
                        tmp = res["derivative"]["smoothed"][0][col]
                        fluor = res["derivative"]["smoothed"][row][col]
                        if np.isfinite(fluor) and fluor > 0.0:
                            bas_temp_min = min(bas_temp_min, float(tmp))
                            bas_temp_max = max(bas_temp_max, float(tmp))
                            bas_fluor_min = min(bas_fluor_min, float(fluor))
                            bas_fluor_max = max(bas_fluor_max, float(fluor))
                            in_bas = [tmp, fluor]
                            bass_json.append(in_bas)
                    # Fixme do not loop over all, use sorted data and clever moving
                    for react in allData["reacts"]:
                        if react["id"] == res["derivative"]["smoothed"][row][0]:
                            for data in react["datas"]:
                                if data["tar"] == res["derivative"]["smoothed"][row][3]:
                                    data["smo"] = list(bass_json)
                allData["smo_temp_min"] = bas_temp_min
                allData["smo_temp_max"] = bas_temp_max
                allData["smo_fluor_min"] = bas_fluor_min
                allData["smo_fluor_max"] = bas_fluor_max

            if "normalized" in res["derivative"]:
                bas_temp_min = 120.0
                bas_temp_max = 0.0
                bas_fluor_min = 99999999
                bas_fluor_max = 0.0
                for row in range(1, len(res["derivative"]["normalized"])):
                    bass_json = []
                    for col in range(6, len(res["derivative"]["normalized"][row])):
                        tmp = res["derivative"]["normalized"][0][col]
                        fluor = res["derivative"]["normalized"][row][col]
                        if np.isfinite(fluor) and fluor > 0.0:
                            bas_temp_min = min(bas_temp_min, float(tmp))
                            bas_temp_max = max(bas_temp_max, float(tmp))
                            bas_fluor_min = min(bas_fluor_min, float(fluor))
                            bas_fluor_max = max(bas_fluor_max, float(fluor))
                            in_bas = [tmp, fluor]
                            bass_json.append(in_bas)
                    # Fixme do not loop over all, use sorted data and clever moving
                    for react in allData["reacts"]:
                        if react["id"] == res["derivative"]["normalized"][row][0]:
                            for data in react["datas"]:
                                if data["tar"] == res["derivative"]["normalized"][row][3]:
                                    data["nrm"] = list(bass_json)
                allData["nrm_temp_min"] = bas_temp_min
                allData["nrm_temp_max"] = bas_temp_max
                allData["nrm_fluor_min"] = bas_fluor_min
                allData["nrm_fluor_max"] = bas_fluor_max

            if "firstDerivative" in res["derivative"]:
                bas_temp_min = 120.0
                bas_temp_max = 0.0
                bas_fluor_min = 99999999
                bas_fluor_max = 0.0
                for row in range(1, len(res["derivative"]["firstDerivative"])):
                    bass_json = []
                    for col in range(6, len(res["derivative"]["firstDerivative"][row])):
                        tmp = res["derivative"]["firstDerivative"][0][col]
                        fluor = res["derivative"]["firstDerivative"][row][col]
                        if np.isfinite(fluor) and fluor > 0.0:
                            bas_temp_min = min(bas_temp_min, float(tmp))
                            bas_temp_max = max(bas_temp_max, float(tmp))
                            bas_fluor_min = min(bas_fluor_min, float(fluor))
                            bas_fluor_max = max(bas_fluor_max, float(fluor))
                            in_bas = [tmp, fluor]
                            bass_json.append(in_bas)
                    # Fixme do not loop over all, use sorted data and clever moving
                    for react in allData["reacts"]:
                        if react["id"] == res["derivative"]["firstDerivative"][row][0]:
                            for data in react["datas"]:
                                if data["tar"] == res["derivative"]["firstDerivative"][row][3]:
                                    data["fdm"] = list(bass_json)
                allData["fdm_temp_min"] = bas_temp_min
                allData["fdm_temp_max"] = bas_temp_max
                allData["fdm_fluor_min"] = bas_fluor_min
                allData["fdm_fluor_max"] = bas_fluor_max

        if "resultsList" in res:
            resList = res["resultsList"]
            for rRow in range(0, len(resList)):
                for rCol in range(0, len(resList[rRow])):
                    if isinstance(resList[rRow][rCol], np.float64) and not np.isfinite(resList[rRow][rCol]):
                        resList[rRow][rCol] = ""
                    if isinstance(resList[rRow][rCol], float) and not math.isfinite(resList[rRow][rCol]):
                        resList[rRow][rCol] = ""
            allData["Meltcurve_Result_Table"] = json.dumps(resList, cls=NpEncoder)

        if "error" in res:
            allData["error"] = res["error"]

        return allData

    def meltCurveAnalysis(self, normMethod="exponential", fluorSource="normalised",
                          truePeakWidth=1.0, artifactPeakWidth=1.0,
                          expoLowTemp=65.0, expoHighTemp=92.0,
                          bilinLowStartTemp=65.0, bilinLowStopTemp=67.0,
                          bilinHighStartTemp=93.0, bilinHighStopTemp=94.0,
                          peakLowTemp=60.0, peakHighTemp=98.0,
                          peakMaxWidth=5.0, peakCutoff=5.0, updateRDML=False,
                          saveRaw=False, saveDerivative=False,
                          saveResultsList=False, saveResultsCSV=False, verbose=False):
        """Performs a melt curve analysis on the run. Modifies the melting temperature values and returns a json with additional data.

        Args:
            self: The class self parameter.
            normMethod: The normalization method "exponential", "bilinear" or "combined"
            fluorSource: choose "normalised" or "smooth" for fluorescence correction factor calculation
            truePeakWidth: the maximum allowed width in Celsius of the expected (true) peak
            artifactPeakWidth: the maximum allowed width in Celsius of all artifact peaks
            expoLowTemp: the low temperature for the exponential normalisation
            expoHighTemp: the high temperature for the exponential normalisation
            bilinLowStartTemp: the low start temperature for the bilinear normalisation
            bilinLowStopTemp: the low stop temperature for the bilinear normalisation
            bilinHighStartTemp: the high start temperature for the bilinear normalisation
            bilinHighStopTemp: the high stop temperature for the bilinear normalisation
            peakLowTemp: peaks below this temperature are ignored
            peakHighTemp: peaks above this temperature are ignored
            peakMaxWidth: peaks broader than this temperature width are ignored
            peakCutoff: the percentage below melting peaks are ignored in calculations
            updateRDML: If true, update the RDML data with the calculated values.
            saveRaw: If true, no raw values are given in the returned data
            saveDerivative: If true, derivative values are given in the returned data
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

        expParent = self._node.getparent()
        rootPar = expParent.getparent()
        dataVersion = rootPar.get('version')
        transSamTar = _sampleTypeToDics(rootPar)

        if dataVersion == "1.0":
            raise RdmlError('MeltCurveAnalysis requires RDML version > 1.0.')
        err = ""

        ##############################
        # Collect the data in arrays #
        ##############################

        # res is a 2 dimensional array accessed only by
        # variables, so columns might be added here
        header = [["id",  # 0
                   "well",  # 1
                   "sample",  # 2
                   "sample type",  # 3
                   "target",   # 4
                   "target chemistry",  # 5
                   "dye",  # 6
                   "excluded",   # 7
                   "note",   # 8
                   "expected melting temperature"]]   # 9
        rar_id = 0
        rar_well = 1
        rar_sample = 2
        rar_sample_type = 3
        rar_tar = 4
        rar_tar_chemistry = 5
        rar_dye = 6
        rar_excl = 7
        rar_note = 8
        rar_exp_melt_temp = 9

        # Never trust user input
        truePeakWidth = float(truePeakWidth)
        artifactPeakWidth = float(artifactPeakWidth)
        expoLowTemp = float(expoLowTemp)
        expoHighTemp = float(expoHighTemp)
        bilinLowStartTemp = float(bilinLowStartTemp)
        bilinLowStopTemp = float(bilinLowStopTemp)
        bilinHighStartTemp = float(bilinHighStartTemp)
        bilinHighStopTemp = float(bilinHighStopTemp)
        peakLowTemp = float(peakLowTemp)
        peakHighTemp = float(peakHighTemp)
        peakMaxWidth = float(peakMaxWidth)
        peakCutoff = float(peakCutoff) / 100.0
        if normMethod not in ["exponential", "bilinear", "combined"]:
            normMethod = "exponential"
        if fluorSource != "normalised":
            fluorSource = "smooth"

        res = []
        finalData = {}
        collAllTemp = {}
        lookUpTemp = {}
        reacts = _get_all_children(self._node, "react")

        # First get the max number of cycles + all temperatures and create the numpy array
        colCount = 0
        for react in reacts:
            react_datas = _get_all_children(react, "data")
            for react_data in react_datas:
                colCount += 1
                mdps = _get_all_children(react_data, "mdp")
                for mdp in mdps:
                    collAllTemp[_get_first_child_text(mdp, "tmp")] = 1

        tempListUnsort = list(collAllTemp.keys())
        tempStrList = sorted(tempListUnsort, key=float)
        tempList = np.array(tempStrList, dtype=np.float64)
        count = 0
        for tp in tempStrList:
            lookUpTemp[tp] = count
            count += 1

        # spFl is the shape for all fluorescence numpy data arrays
        spFl = (colCount, len(tempList))
        rawFluor = np.zeros(spFl, dtype=np.float64)
        rawFluor[rawFluor < 1] = np.nan

        # Initialization of the vecNoAmplification vector
        vecExcludedByUser = np.zeros(spFl[0], dtype=np.bool_)
        rdmlElemData = []

        # Now process the data for numpy and create results array
        anyRawData = False
        rowCount = 0
        for react in reacts:
            posId = react.get('id')
            pIdNumber = (int(posId) - 1) % int(self["pcrFormat_columns"]) + 1
            pIdLetter = chr(ord("A") + int((int(posId) - 1) / int(self["pcrFormat_columns"])))
            pWell = pIdLetter + str(pIdNumber)
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
                excl = _get_first_child_text(react_data, "excl")
                if not excl == "":
                    vecExcludedByUser[rowCount] = True
                noteVal = _get_first_child_text(react_data, "note")
                rdmlElemData.append(react_data)
                res.append([posId, pWell, sample, "",  target, "", "",  excl, noteVal, ""])  # Must match header length
                mdps = _get_all_children(react_data, "mdp")
                for mdp in mdps:
                    cTemp = _get_first_child_text(mdp, "tmp")
                    cFluor = _get_first_child_text(mdp, "fluor")
                    if cTemp != "" and cFluor != "" and cTemp in lookUpTemp:
                        rawFluor[rowCount, lookUpTemp[cTemp]] = float(cFluor)
                        anyRawData = True
                rowCount += 1
        if anyRawData == False:
            raise RdmlError('Melting Curve Analysis requires raw data. No raw data were found in this run.')

        # Look up sample and target information
        parExp = self._node.getparent()
        parRoot = parExp.getparent()

        dicLU_dyes = {}
        luDyes = _get_all_children(parRoot, "dye")
        for lu_dye in luDyes:
            lu_chemistry = _get_first_child_text(lu_dye, "dyeChemistry")
            if lu_chemistry == "":
                lu_chemistry = "non-saturating DNA binding dye"
            if lu_dye.attrib['id'] != "":
                dicLU_dyes[lu_dye.attrib['id']] = lu_chemistry

        dicLU_targets = {}
        dicLU_tarMelt = {}
        luTargets = _get_all_children(parRoot, "target")
        for lu_target in luTargets:
            meltingTemperature = _get_first_child_text(lu_target, "meltingTemperature")
            forId = _get_first_child(lu_target, "dyeId")
            lu_dyeId = ""
            if forId is not None:
                if forId.attrib['id'] != "":
                    lu_dyeId = forId.attrib['id']
            if lu_dyeId == "" or lu_dyeId not in dicLU_dyes:
                dicLU_targets[lu_target.attrib['id']] = "non-saturating DNA binding dye"
            if lu_target.attrib['id'] != "":
                dicLU_targets[lu_target.attrib['id']] = dicLU_dyes[lu_dyeId]
                dicLU_tarMelt[lu_target.attrib['id']] = meltingTemperature

        # Update the table with dictionary help
        for oRow in range(0, spFl[0]):
            if res[oRow][rar_sample] != "":
                if res[oRow][rar_tar] != "":
                    res[oRow][rar_sample_type] = transSamTar[res[oRow][rar_sample]][res[oRow][rar_tar]]
            if res[oRow][rar_tar] != "":

                res[oRow][rar_tar_chemistry] = dicLU_targets[res[oRow][rar_tar]]
                res[oRow][rar_exp_melt_temp] = dicLU_tarMelt[res[oRow][rar_tar]]

        if saveRaw:
            rawData = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar],
                        header[0][rar_excl], header[0][rar_exp_melt_temp]]]
            for oCol in tempStrList:
                rawData[0].append(oCol)
            for oRow in range(0, spFl[0]):
                rawData.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar],
                                res[oRow][rar_excl], res[oRow][rar_exp_melt_temp]])
                for oCol in range(0, spFl[1]):
                    rawData[oRow + 1].append(float(rawFluor[oRow, oCol]))
            finalData["rawData"] = rawData

        # Count the targets and create the target variables
        # Position 0 is for the general over all window without targets
        vecTarget = np.zeros(spFl[0], dtype=np.int64)
        vecTarget[vecTarget <= 0] = -1
        targetsCount = 1
        tarWinLookup = {}
        tarReverseLookup = {}
        for oRow in range(0, spFl[0]):
            if res[oRow][rar_tar] not in tarWinLookup:
                tarWinLookup[res[oRow][rar_tar]] = targetsCount
                tarReverseLookup[targetsCount] = res[oRow][rar_tar]
                targetsCount += 1
            vecTarget[oRow] = tarWinLookup[res[oRow][rar_tar]]
        bilinLowStart = np.zeros(targetsCount, dtype=np.float64)
        bilinLowStop = np.zeros(targetsCount, dtype=np.float64)
        bilinHighStart = np.zeros(targetsCount, dtype=np.float64)
        bilinHighStop = np.zeros(targetsCount, dtype=np.float64)

        LowTm = np.zeros(targetsCount, dtype=np.float64)
        HighTm = np.zeros(targetsCount, dtype=np.float64)

        # Initialization of the error vectors
        vecSweeps = np.zeros(spFl[0], dtype=np.bool_)

        #########################
        # Get the data in shape #
        #########################

        # There should be no negative values in uncorrected raw data
        absMinFluor = np.nanmin(rawFluor)
        if absMinFluor < 0.0:
            err = _add_err(err, ";", "Error: Fluorescence data have negative values. Use raw data without baseline correction! ")
            err = _add_err(err, "", "Baseline corrected data not using a constant factor will result in wrong melting curves!")

        # Initial smooth of raw data
        smoothFluor = _mca_smooth(tempList, rawFluor)

        # Exponential normalisation
        normalMelting = np.zeros(spFl, dtype=np.float64)
        posLowT = 0
        while posLowT < spFl[1] - 1 and tempList[posLowT] < expoLowTemp:
            posLowT += 1
        posHighT = spFl[1] - 1
        while posHighT > 0 and tempList[posHighT] > expoHighTemp:
            posHighT -= 1

        for rRow in range(0, spFl[0]):  # loop rRow for every reaction
            # calculate FDLow and FDHigh from mc[]
            FDLow = -1 * (smoothFluor[rRow][posLowT] - smoothFluor[rRow][posLowT - 1])
            FDHigh = -1 * (smoothFluor[rRow][posHighT] - smoothFluor[rRow][posHighT - 1])

            # Rarely happens, protects the log from negative values
            if FDLow <= 0.0:
                FDLow = 0.00001
            if FDHigh <= 0.0:
                FDHigh = 0.000001  # not same as FDLow = 0.00001

            # determine Aexp and Cexp
            Aexp = (np.log(FDHigh) - np.log(FDLow)) / (expoHighTemp - expoLowTemp)
            Cexp = -1 * FDLow / Aexp

            # apply exponential base trend correction
            MaxMCCorr = 0.0
            MinMCCorr = 10000.0
            for rCol in range(0, spFl[1]):
                normalMelting[rRow][rCol] = smoothFluor[rRow][rCol] - Cexp * np.exp(Aexp * (tempList[rCol] - expoLowTemp))
                if normalMelting[rRow][rCol] > MaxMCCorr:
                    MaxMCCorr = normalMelting[rRow][rCol]
                if normalMelting[rRow][rCol] < MinMCCorr:
                    MinMCCorr = normalMelting[rRow][rCol]

            for rCol in range(0, spFl[1]):
                normalMelting[rRow][rCol] = (normalMelting[rRow][rCol] - MinMCCorr) / (MaxMCCorr - MinMCCorr)

        if normMethod in ["bilinear", "combined"]:
            ##################################
            # Finding the suitable low range #
            ##################################
            # determine index of low start temperature
            startindex = 0
            while tempList[startindex] < bilinLowStartTemp and startindex < len(tempList) - 1:
                startindex += 1
            # determine index of low stop temperature
            stopindex = len(tempList) - 1
            while tempList[stopindex] > bilinLowStopTemp and stopindex > 0:
                stopindex -= 1

            NtempsInRange = stopindex - startindex + 1

            # determine index of high start temperature
            starthighT = 0
            while tempList[starthighT] < bilinHighStartTemp and starthighT < len(tempList) - 1:
                starthighT += 1
            # determine index of high stop temperature
            stophighT = len(tempList) - 1
            while tempList[stophighT] > bilinHighStopTemp and stophighT > 0:
                stophighT -= 1

            MeanSlope = np.zeros((targetsCount, 3 * NtempsInRange + 1), dtype=np.float64)
            SDSlope = np.zeros((targetsCount, 3 * NtempsInRange + 1), dtype=np.float64)
            IndexR = -1
            for k in range(startindex, startindex + 3 * NtempsInRange + 1):
                startlowT = k - 1
                stoplowT = startlowT + NtempsInRange
                IndexR += 1  # counts number of tested T ranges
                SumSlopes = np.zeros(targetsCount, dtype=np.float64)
                SumSlopes2 = np.zeros(targetsCount, dtype=np.float64)
                cntSlopes = np.zeros(targetsCount, dtype=np.int64)

                if normMethod == "combined":
                    bilinNormal = normalMelting.copy()
                else:
                    bilinNormal = smoothFluor.copy()

                [slopelow, interceptlow] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, startlowT, stoplowT)
                [slopehigh, intercepthigh] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, starthighT, stophighT)
                LowTline = interceptlow[:, np.newaxis] + slopelow[:, np.newaxis] * tempList
                HighTline = intercepthigh[:, np.newaxis] + slopehigh[:, np.newaxis] * tempList
                bilinNormal = (bilinNormal - HighTline) / (LowTline - HighTline)
                [slopeNMC, interceptNMC] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal,
                                                       startlowT + NtempsInRange, startlowT + 2 * NtempsInRange)

                for j in range(0, spFl[0]):
                    nonSweep = True
                    for i in range(0, spFl[1] - 1):
                        if (np.abs(bilinNormal[j][i] - bilinNormal[j][i + 1]) > 0.1 and
                                (bilinNormal[j][i] > 1.0 > bilinNormal[j][i + 1] or
                                 bilinNormal[j][i] < 1.0 < bilinNormal[j][i + 1])):
                            nonSweep = False
                    if nonSweep:
                        curTarNr = tarWinLookup[res[j][rar_tar]]
                        SumSlopes[curTarNr] += slopeNMC[j]
                        SumSlopes2[curTarNr] += slopeNMC[j] * slopeNMC[j]
                        cntSlopes[curTarNr] += 1

                for curTarNr in range(1, targetsCount):
                    if cntSlopes[curTarNr] > 1:
                        MeanSlope[curTarNr][IndexR] = SumSlopes[curTarNr] / cntSlopes[curTarNr]
                        SDSlope[curTarNr][IndexR] = np.sqrt((SumSlopes2[curTarNr] - (SumSlopes[curTarNr] *
                                                            SumSlopes[curTarNr] / cntSlopes[curTarNr])) /
                                                            (cntSlopes[curTarNr] - 1))
                    else:
                        if IndexR == 0:
                            MeanSlope[curTarNr][IndexR] = 10.0
                            SDSlope[curTarNr][IndexR] = 0.0
                        else:
                            MeanSlope[curTarNr][IndexR] = MeanSlope[curTarNr][IndexR - 1]
                            SDSlope[curTarNr][IndexR] = SDSlope[curTarNr][IndexR - 1]

            MinSlope = 10.0  # default set to 10.0
            IndexMin = 0

            for curTarNr in range(1, targetsCount):
                for k in range(0, IndexR - 1):
                    if normMethod == "combined":
                        if SDSlope[curTarNr][k] < MinSlope:
                            MinSlope = SDSlope[curTarNr][k]
                            IndexMin = k
                    else:
                        CritSlope = MeanSlope[curTarNr][k] + 2 * SDSlope[curTarNr][k]
                        if CritSlope < 0.0 and (0.0 - CritSlope) < MinSlope:
                            MinSlope = 0.0 - CritSlope
                            IndexMin = k
                bilinLowStart[curTarNr] = tempList[startindex + IndexMin]
                bilinLowStop[curTarNr] = tempList[startindex + IndexMin + NtempsInRange]
                LowTm[curTarNr] = bilinLowStop[curTarNr]

            ###################################
            # Finding the suitable high range #
            ###################################
            for curTarNr in range(1, targetsCount):
                starttemp = bilinLowStart[curTarNr]
                # determine index of low start temperature
                startlowT = 0
                while tempList[startlowT] < starttemp and startlowT < len(tempList) - 1:
                    startlowT += 1

                stoptemp = bilinLowStop[curTarNr]
                # determine index of low stop temperature
                stoplowT = len(tempList) - 1
                while tempList[stoplowT] > stoptemp and stoplowT > 0:
                    stoplowT -= 1

                # determine startindex for range of high temps to test
                starttemp = 90.0
                # determine indices of high start and stop temperature
                startindex = 1
                while tempList[startindex] < starttemp and startindex < len(tempList) - 1:
                    startindex += 1

                stoptemp = 91.0
                # determine indices of high start and stop temperature
                stopindex = len(tempList) - 1
                while tempList[stopindex] > stoptemp and stopindex > 0:
                    stopindex -= 1

                NtempsInRange = stopindex - startindex + 1

                MeanVal = np.zeros((targetsCount, 5 * NtempsInRange + 1), dtype=np.float64)
                SDVal = np.zeros((targetsCount, 5 * NtempsInRange + 1), dtype=np.float64)
                IndexR = -1
                for k in range(startindex, startindex + 5 * NtempsInRange + 1):
                    starthighT = k - 1
                    stophighT = starthighT + NtempsInRange
                    IndexR += 1  # counts number of tested T ranges

                    SumVal = 0.0
                    SumVal2 = 0.0
                    cntVal = 0

                    if normMethod == "combined":
                        bilinNormal = normalMelting.copy()
                    else:
                        bilinNormal = smoothFluor.copy()

                    [slopelow, interceptlow] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, startlowT, stoplowT)
                    [slopehigh, intercepthigh] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, starthighT, stophighT)
                    LowTline = interceptlow[:, np.newaxis] + slopelow[:, np.newaxis] * tempList
                    HighTline = intercepthigh[:, np.newaxis] + slopehigh[:, np.newaxis] * tempList
                    bilinNormal = (bilinNormal - HighTline) / (LowTline - HighTline)
                    # [slopeNMC, interceptNMC] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal,
                    #                                        startlowT + NtempsInRange, startlowT + 2 * NtempsInRange)

                    for j in range(0, spFl[0]):
                        if curTarNr == tarWinLookup[res[j][rar_tar]]:
                            nonSweep = True
                            for i in range(k - 2 * NtempsInRange, k + 1):
                                if (np.abs(bilinNormal[j][i] - bilinNormal[j][i + 1]) > 0.1 and
                                        (bilinNormal[j][i] > 0.0 > bilinNormal[j][i + 1] or
                                         bilinNormal[j][i] < 0.0 < bilinNormal[j][i + 1])):
                                    nonSweep = False
                            if nonSweep:
                                sumValIn = 0.0
                                cntValIn = 0
                                for i in range(k - NtempsInRange, k + 1):
                                    sumValIn += bilinNormal[j][i]
                                    cntValIn += 1
                                MeanValIn = sumValIn / cntValIn
                                SumVal += MeanValIn
                                SumVal2 += MeanValIn * MeanValIn
                                cntVal += 1

                    MeanVal[curTarNr][IndexR] = SumVal / cntVal
                    SDVal[curTarNr][IndexR] = np.sqrt((SumVal2 - (SumVal * SumVal / cntVal)) / (cntVal - 1))

                MaxVal = -10.0  # default set to -10.0
                IndexMax = 0

                # criterion is the same for bi-linear or exponential+bi-linear
                for k in range(0, IndexR):
                    CritVal = MeanVal[curTarNr][k] - 2 * SDVal[curTarNr][k]
                    if CritVal > MaxVal:
                        MaxVal = CritVal
                        IndexMax = k

                bilinHighStart[curTarNr] = tempList[startindex + IndexMax]
                bilinHighStop[curTarNr] = tempList[startindex + IndexMax + NtempsInRange]
                HighTm[curTarNr] = bilinHighStart[curTarNr]

            #################################
            # Do the bilinear normalisation #
            #################################
            if normMethod == "combined":
                bilinNormal = normalMelting.copy()
            else:
                bilinNormal = smoothFluor.copy()

            for curTarNr in range(1, targetsCount):
                # determine index of low start temperature
                startlowT = 0
                while tempList[startlowT] < bilinLowStart[curTarNr] and startlowT < len(tempList) - 1:
                    startlowT += 1
                # determine index of low stop temperature
                stoplowT = len(tempList) - 1
                while tempList[stoplowT] > bilinLowStop[curTarNr] and stoplowT > 0:
                    stoplowT -= 1
                # determine indices of high start and stop temperature
                starthighT = 1
                while tempList[starthighT] < bilinHighStart[curTarNr] and starthighT < len(tempList) - 1:
                    starthighT += 1
                # determine indices of high start and stop temperature
                stophighT = len(tempList) - 1
                while tempList[stophighT] > bilinHighStop[curTarNr] and stophighT > 0:
                    stophighT -= 1

                [slopelow, interceptlow] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, startlowT, stoplowT)
                [slopehigh, intercepthigh] = _mca_linReg(np.tile(tempList, (spFl[0], 1)), bilinNormal, starthighT, stophighT)

                LowTline = interceptlow[:, np.newaxis] + slopelow[:, np.newaxis] * tempList
                HighTline = intercepthigh[:, np.newaxis] + slopehigh[:, np.newaxis] * tempList

                for j in range(0, spFl[0]):
                    if curTarNr == tarWinLookup[res[j][rar_tar]]:
                        normalMelting[j] = (bilinNormal[j] - HighTline[j]) / (LowTline[j] - HighTline[j])
                        for i in range(0, spFl[1]):
                            # avoid sweeps because LowTline and HighTline are about to cross
                            if i > stophighT:
                                if abs(normalMelting[j][i] - normalMelting[j][i - 1]) > 1.01 * normalMelting[j][i - 1]:
                                    normalMelting[j][i] = normalMelting[j][i - 1]

        # FindSweepsButtonClick ???

        # Derivate normalMelting
        tmp = (tempList + np.roll(tempList, 1)) / 2  # Shift to right
        rawFirstDerivativeTemp = tmp[1:]
        tmp = np.roll(normalMelting, 1, axis=1) - normalMelting  # Shift to right
        rawFirstDerivative = tmp[:, 1:]

        # Delete the first three columns
        rawFirstDerivativeTemp = rawFirstDerivativeTemp[3:]
        rawFirstDerivative = rawFirstDerivative[:, 3:]

        # Smooth of raw data
        smoothFirstDerivative = _mca_smooth(rawFirstDerivativeTemp, rawFirstDerivative)

        # Derivate smoothFirstDerivative
        rawSecondDerivativeTemp = tempList[4:-2]
        tmp = smoothFirstDerivative - np.roll(smoothFirstDerivative, 1, axis=1)  # Shift to right
        rawSecondDerivative = tmp[:, 1:-1]

        # Smooth of raw data
        smoothSecondDerivative = _mca_smooth(rawSecondDerivativeTemp, rawSecondDerivative)

        # Save derivative data
        if saveDerivative:
            finalData["derivative"] = {}

            smoothData = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar],
                         header[0][rar_excl], header[0][rar_exp_melt_temp]]]
            for oCol in tempStrList:
                smoothData[0].append(oCol)
            for oRow in range(0, spFl[0]):
                smoothData.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar],
                                res[oRow][rar_excl], res[oRow][rar_exp_melt_temp]])
                for oCol in range(0, spFl[1]):
                    smoothData[oRow + 1].append(float(smoothFluor[oRow, oCol]))
            finalData["derivative"]["smoothed"] = smoothData

            normData = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar],
                         header[0][rar_excl], header[0][rar_exp_melt_temp]]]
            for oCol in tempStrList:
                normData[0].append(oCol)
            for oRow in range(0, spFl[0]):
                normData.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar],
                                 res[oRow][rar_excl], res[oRow][rar_exp_melt_temp]])
                for oCol in range(0, spFl[1]):
                    normData[oRow + 1].append(float(normalMelting[oRow, oCol]))
            finalData["derivative"]["normalized"] = normData

            firstDerData = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar],
                         header[0][rar_excl], header[0][rar_exp_melt_temp]]]
            for oCol in rawFirstDerivativeTemp:
                firstDerData[0].append(oCol)
            for oRow in range(0, spFl[0]):
                firstDerData.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar],
                                     res[oRow][rar_excl], res[oRow][rar_exp_melt_temp]])
                for oCol in range(0, smoothFirstDerivative.shape[1]):
                    firstDerData[oRow + 1].append(float(smoothFirstDerivative[oRow, oCol]))
            finalData["derivative"]["firstDerivative"] = firstDerData

            secondDerData = [[header[0][rar_id], header[0][rar_well], header[0][rar_sample], header[0][rar_tar],
                              header[0][rar_excl], header[0][rar_exp_melt_temp]]]
            for oCol in rawSecondDerivativeTemp:
                secondDerData[0].append(oCol)
            for oRow in range(0, spFl[0]):
                secondDerData.append([res[oRow][rar_id], res[oRow][rar_well], res[oRow][rar_sample], res[oRow][rar_tar],
                                res[oRow][rar_excl], res[oRow][rar_exp_melt_temp]])
                for oCol in range(0, smoothSecondDerivative.shape[1]):
                    secondDerData[oRow + 1].append(float(smoothSecondDerivative[oRow, oCol]))
            finalData["derivative"]["secondDerivative"] = secondDerData

        #######################################
        # Now find peaks and their parameters #
        #######################################
        if saveResultsList:
            peakResTemp = []
            peakResWidth = []
            peakResH = []
            peakResSumH = []
            peakResDeltaH = []
            peakResSumDeltaH = []
            peakResFluor = []
            peakResSumFuor = []
            truePeakFinPos = [-1] * spFl[0]
            for pos in range(0, spFl[0]):  # loop rRow for every reaction
                peakResTemp.append([])
                peakResWidth.append([])
                peakResH.append([])
                peakResSumH.append(0.0)
                peakResDeltaH.append([])
                peakResSumDeltaH.append(0.0)
                peakResFluor.append([])
                peakResSumFuor.append(0.0)
                for fdPos in range(1, len(rawFirstDerivativeTemp) - 2):
                    if smoothFirstDerivative[pos][fdPos - 1] <= smoothFirstDerivative[pos][fdPos] > smoothFirstDerivative[pos][fdPos + 1]:
                        # found a peak
                        peakTemp = rawFirstDerivativeTemp[fdPos]
                        sdPosMax = fdPos - 1
                        sdValMax = smoothSecondDerivative[pos][sdPosMax]
                        sdPosMin = fdPos
                        sdValMin = smoothSecondDerivative[pos][sdPosMin]
                        while sdPosMax > 0 and sdValMax < smoothSecondDerivative[pos][sdPosMax - 1]:
                            sdPosMax -= 1
                            sdValMax = smoothSecondDerivative[pos][sdPosMax]
                        while (sdPosMin < (len(rawSecondDerivativeTemp) - 2) and
                               sdValMin > smoothSecondDerivative[pos][sdPosMin + 1]):
                            sdPosMin += 1
                            sdValMin = smoothSecondDerivative[pos][sdPosMin]

                        if sdPosMax > 0:  # All peaks have to go up first!
                            lowPeakTemp = rawSecondDerivativeTemp[sdPosMax]

                            ndLowPeakPos = sdPosMax + 4  # is the precise same temp
                            if fluorSource == "normalised":
                                lowFinFluor = normalMelting[pos][ndLowPeakPos]
                            else:
                                lowFinFluor = smoothFluor[pos][ndLowPeakPos]

                            validPeakSD = True
                            if sdPosMin > (len(rawSecondDerivativeTemp) - 4):
                                # assume symmetry when down side is missing
                                peakTempWitdth = 2 * (peakTemp - lowPeakTemp)
                                ndPeakPos = fdPos + 3  # The peak is between +3 and +4
                                if fluorSource == "normalised":
                                    medFinFluor = normalMelting[pos][ndPeakPos]
                                else:
                                    medFinFluor = smoothFluor[pos][ndPeakPos]
                                fluorDrop = 2 * (lowFinFluor - medFinFluor)

                                # sdPosMax position is on FD between sdPosMax and sdPosMax + 1
                                deltaH = smoothFirstDerivative[pos][fdPos] - smoothFirstDerivative[pos][sdPosMax]
                                if (smoothFirstDerivative[pos][fdPos] < smoothFirstDerivative[pos][sdPosMax] or
                                    smoothFirstDerivative[pos][fdPos] < 0.0 or
                                    smoothFirstDerivative[pos][sdPosMax] < 0.0):
                                    validPeakSD = False
                            else:
                                highPeakTemp = rawSecondDerivativeTemp[sdPosMin]
                                peakTempWitdth = highPeakTemp - lowPeakTemp
                                ndHighPeakPos = sdPosMin + 4  # is the precise same temp
                                if fluorSource == "normalised":
                                    highFinFluor = normalMelting[pos][ndHighPeakPos]
                                else:
                                    highFinFluor = smoothFluor[pos][ndHighPeakPos]
                                fluorDrop = lowFinFluor - highFinFluor

                                # sdPosMax position is on FD between sdPosMax and sdPosMax + 1
                                if smoothFirstDerivative[pos][sdPosMin] < 0.0:
                                    deltaH = smoothFirstDerivative[pos][fdPos] - smoothFirstDerivative[pos][sdPosMax]
                                else:
                                    lowH = smoothFirstDerivative[pos][sdPosMax] + smoothFirstDerivative[pos][sdPosMin + 1]
                                    deltaH = smoothFirstDerivative[pos][fdPos] - lowH / 2
                                if (smoothFirstDerivative[pos][fdPos] < smoothFirstDerivative[pos][sdPosMax] or
                                    smoothFirstDerivative[pos][fdPos] < smoothFirstDerivative[pos][sdPosMin + 1] or
                                    smoothFirstDerivative[pos][fdPos] < 0.0 or
                                    smoothFirstDerivative[pos][sdPosMax] < 0.0):
                                    validPeakSD = False

                            if (fluorDrop > 0.0 and
                                    peakTempWitdth < peakMaxWidth and
                                    peakLowTemp < peakTemp < peakHighTemp and
                                    validPeakSD):
                                peakResTemp[pos].append(peakTemp)
                                peakResWidth[pos].append(peakTempWitdth)
                                peakResH[pos].append(smoothFirstDerivative[pos][fdPos])
                                peakResSumH[pos] += smoothFirstDerivative[pos][fdPos]
                                peakResDeltaH[pos].append(deltaH)
                                peakResSumDeltaH[pos] += deltaH
                                peakResFluor[pos].append(fluorDrop)
                                peakResSumFuor[pos] += fluorDrop

            # Set unwanted peaks below peakCutoff to -10.0
            for oRow in range(0, spFl[0]):
                for oCol in range(0, len(peakResTemp[oRow])):
                    if peakCutoff > peakResDeltaH[oRow][oCol] / peakResSumDeltaH[oRow]:
                        peakResTemp[oRow][oCol] = -10.0

            # Recalculate the sums
            for oRow in range(0, spFl[0]):
                peakResSumH[oRow] = 0.0
                peakResSumDeltaH[oRow] = 0.0
                peakResSumFuor[oRow] = 0.0
                for oCol in range(0, len(peakResTemp[oRow])):
                    if peakResTemp[oRow][oCol] > 0.0:
                        peakResSumH[oRow] += peakResH[oRow][oCol]
                        peakResSumDeltaH[oRow] += peakResDeltaH[oRow][oCol]
                        peakResSumFuor[oRow] += peakResFluor[oRow][oCol]

            # Find the expected peak
            checkedPeakTemp = [row[:] for row in peakResTemp]
            expTemp = [[]]
            for curTarNr in range(1, targetsCount):
                expTemp.append(dicLU_tarMelt[tarReverseLookup[curTarNr]])
                if expTemp[curTarNr] == "" or float(expTemp[curTarNr]) < 20.0 or float(expTemp[curTarNr]) > 100.0:
                    expTemp[curTarNr] = -10.0
                    # if False:
                    #     # Find the best peak as true peak
                    #     startPeakTemp = tempList[0] + truePeakWidth
                    #     stopPeakTemp = tempList[-1] - truePeakWidth
                    #     startPeakPos = 0
                    #     stopPeakPos = len(tempList) - 2
                    #     while tempList[startPeakPos + 1] < startPeakTemp:
                    #         startPeakPos += 1
                    #     while tempList[stopPeakPos - 1] > stopPeakTemp:
                    #         stopPeakPos -= 1
                    #     maxDeltaH = 0.0
                    #     for curTempPos in range(startPeakPos, stopPeakPos + 1):
                    #         curTemp = tempList[curTempPos]
                    #         curPeakInRange = 0
                    #         meanCalc = 0.0
                    #         deltaHSum = 0.0
                    #         for oRow in range(0, spFl[0]):
                    #             if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                    #                 for oCol in range(0, len(checkedPeakTemp[oRow])):
                    #                     if curTemp - truePeakWidth < checkedPeakTemp[oRow][oCol] < curTemp + truePeakWidth:
                    #                         curPeakInRange += 1
                    #                         meanCalc += checkedPeakTemp[oRow][oCol]
                    #                         deltaHSum += peakResDeltaH[oRow][oCol]
                    #         # if curPeakInRange >= maxPeakInRange and curPeakInRange > 0.0:
                    #         if curPeakInRange > 0.0 and maxDeltaH < deltaHSum / curPeakInRange:
                    #             maxDeltaH = deltaHSum / curPeakInRange
                    #             expTemp[curTarNr] = meanCalc / curPeakInRange
                # print(str(tarReverseLookup[curTarNr]) + " - " + str(expTemp[curTarNr]))

                # Now we have a peak temp
                if float(expTemp[curTarNr]) > 0.0:
                    for oRow in range(0, spFl[0]):
                        if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                            maxDeltaH = 0.0
                            for oCol in range(0, len(checkedPeakTemp[oRow])):
                                if float(expTemp[curTarNr]) - truePeakWidth < checkedPeakTemp[oRow][oCol] < float(expTemp[curTarNr]) + truePeakWidth:
                                    if peakResDeltaH[oRow][oCol] >= maxDeltaH:
                                        # Choose the peak with the biggest deltaH
                                        maxDeltaH = peakResDeltaH[oRow][oCol]
                                        truePeakFinPos[oRow] = oCol
                            if maxDeltaH > 0.0:
                                checkedPeakTemp[oRow][truePeakFinPos[oRow]] = -10.0

            artifactPeaks = [[]]
            orderPeaks = -1 * np.ones((spFl[0], 1), dtype=np.int64)
            for curTarNr in range(1, targetsCount):
                # Find the artifact peaks
                startPeakTemp = tempList[0] + artifactPeakWidth
                stopPeakTemp = tempList[-1] - artifactPeakWidth
                startPeakPos = 0
                stopPeakPos = len(tempList) - 2
                while tempList[startPeakPos + 1] < startPeakTemp:
                    startPeakPos += 1
                while tempList[stopPeakPos - 1] > stopPeakTemp:
                    stopPeakPos -= 1
                availPos = list(range(startPeakPos, stopPeakPos + 1))

                # Remove target temp from list
                delPos = -1
                for listPos in range(0, len(availPos)):
                    if tempList[availPos[listPos]] == expTemp[curTarNr]:
                        delPos = listPos
                if delPos > -1:
                    del availPos[delPos]

                artifactPeaks.append([])
                stillPeaks = True
                peakCount = -1
                while stillPeaks:
                    stillPeaks = False
                    foundPeaks = False
                    maxPeakInRange = 0
                    maxDeltaHSum = 0.0
                    goodTempPos = [-1.0]
                    goodPosList = -1 * np.ones((spFl[0], 1), dtype=np.int64)
                    goodPosCount = 0
                    for curTempPos in availPos:
                        curTemp = tempList[curTempPos]
                        curPeakInRange = 0
                        curDeltaHSum = 0.0
                        curPosList = -1 * np.ones((spFl[0], 1), dtype=np.int64)
                        for oRow in range(0, spFl[0]):
                            if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                                maxDeltaH = 0.0
                                curGoodCol = -1
                                for oCol in range(0, len(checkedPeakTemp[oRow])):
                                    if curTemp - artifactPeakWidth < checkedPeakTemp[oRow][oCol] < curTemp + artifactPeakWidth:
                                        if peakResDeltaH[oRow][oCol] >= maxDeltaH:
                                            # Choose the peak with the biggest deltaH
                                            maxDeltaH = peakResDeltaH[oRow][oCol]
                                            curGoodCol = oCol
                                if maxDeltaH > 0.0:
                                    foundPeaks = True
                                    curPeakInRange += 1
                                    curDeltaHSum += maxDeltaH
                                    curPosList[oRow] = curGoodCol

                        if curPeakInRange >= maxPeakInRange and curPeakInRange > 0.0:
                            if curDeltaHSum > maxDeltaHSum or curPeakInRange > maxPeakInRange:
                                goodTempPos = [curTempPos]
                                goodPosList = curPosList
                                goodPosCount = 0
                                maxDeltaHSum = curDeltaHSum
                                maxPeakInRange = curPeakInRange
                            else:
                                if curDeltaHSum == maxDeltaHSum or curPeakInRange == maxPeakInRange:
                                    goodTempPos.append(curTempPos)
                                    goodPosList = np.concatenate((goodPosList, curPosList), axis=1)
                                    goodPosCount += 1

                    if foundPeaks and maxDeltaHSum > 0.0:
                        peakCount += 1  # On pos 0 is expected peak
                        if peakCount >= len(orderPeaks[0]):
                            ttPeaks = -1 * np.ones((spFl[0], 1), dtype=np.int64)
                            orderPeaks = np.concatenate((orderPeaks, ttPeaks), axis=1)

                        usePos = int(goodPosCount / 2)
                        stillPeaks = True
                        artifactPeaks[curTarNr].append(tempList[goodTempPos[usePos]])
                        # Avoid duplicate temps
                        del availPos[availPos.index(goodTempPos[usePos])]
                        for oRow in range(0, spFl[0]):
                            if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                                colPos = goodPosList[oRow][usePos]
                                if colPos >= 0:
                                    orderPeaks[oRow][peakCount] = colPos
                                    checkedPeakTemp[oRow][colPos] = -10.0
                                    # add to result array

            # for oRow in range(0, len(checkedPeakTemp)):
            #     for oCol in range(0, len(checkedPeakTemp[oRow])):
            #         if checkedPeakTemp[oRow][colPos] > -10.0:
            #             print("Error: Peak not reported in Row " + str(oRow) + " Col " + str(oCol) + " with: " +
            #                   str(checkedPeakTemp[oRow][colPos]))

            # Calculate the dimensions
            calcRows = 1 + (targetsCount - 1) + spFl[0]  # header, sum cols, data Cols
            maxPeaks = 0
            for curTarNr in range(1, targetsCount):
                peakCount = len(artifactPeaks[curTarNr])
                if float(expTemp[curTarNr]) > 0.0:
                    peakCount += 1
                if peakCount > maxPeaks:
                    maxPeaks = peakCount
            calcCols = 8 + 1 + 8 + 1 + maxPeaks * (8 + 1) - 1  # rdml data, good peak, all peaks
            # Create the table
            resTable = [["" for x in range(calcCols)] for y in range(calcRows)]

            resTable[0][0] = header[0][rar_id]
            resTable[0][1] = header[0][rar_well]
            resTable[0][2] = header[0][rar_sample]
            resTable[0][3] = header[0][rar_sample_type]
            resTable[0][4] = header[0][rar_tar]
            resTable[0][5] = header[0][rar_excl]
            resTable[0][6] = header[0][rar_note]
            resTable[0][7] = header[0][rar_exp_melt_temp]

            for count in range(0, maxPeaks + 1):
                offset = 8 + count * 9
                resTable[0][offset] = "peak temp"
                resTable[0][offset + 1] = "peak width"
                resTable[0][offset + 2] = "Fluor"
                resTable[0][offset + 3] = "peak correction factor"
                resTable[0][offset + 4] = "deltaH"
                resTable[0][offset + 5] = "deltaH factor"
                resTable[0][offset + 6] = "H"
                resTable[0][offset + 7] = "H factor"

            rowPos = 0
            for curTarNr in range(1, targetsCount):
                # Prepare average row
                rowPos += 1
                averageRow = rowPos
                resTable[averageRow][0] = "--"
                resTable[averageRow][1] = 0
                resTable[averageRow][2] = "--"
                resTable[averageRow][3] = "Average"
                resTable[averageRow][4] = tarReverseLookup[curTarNr]
                resTable[averageRow][5] = "--"
                resTable[averageRow][6] = "--"
                resTable[averageRow][7] = dicLU_tarMelt[tarReverseLookup[curTarNr]]

                # Create the true peak columns
                meanResTemp = 0.0
                meanResPeakWidth = 0.0
                meanResH = 0.0
                meanResSumH = 0.0
                meanResDeltaH = 0.0
                meanResSumDeltaH = 0.0
                meanResFluor = 0.0
                meanResSumFuor = 0.0
                countAddedRows = 0
                for oRow in range(0, spFl[0]):
                    if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                        rowPos += 1
                        ###################################
                        # calculate excl and note strings #
                        ###################################
                        exclVal = _cleanErrorString(res[oRow][rar_excl], "melt")
                        noteVal = _cleanErrorString(res[oRow][rar_note], "melt")
                        finPeakCount = 0

                        for curPeak in peakResTemp[oRow]:
                            if curPeak > 0.0:
                                finPeakCount += 1

                        if res[oRow][rar_exp_melt_temp] != "":
                            if truePeakFinPos[oRow] < 0.0:
                                if res[oRow][rar_sample_type] in ["std", "pos"]:
                                    exclVal += "no product with expected melting temperature;"
                                if res[oRow][rar_sample_type] in ["unkn"]:
                                    noteVal += "no product with expected melting temperature;"
                                if finPeakCount == 1:
                                    if res[oRow][rar_sample_type] in ["ntc", "nac", "ntp", "nrt"]:
                                        noteVal += "product with different melting temperatures detected;"
                                    if res[oRow][rar_sample_type] in ["std", "pos", "unkn"]:
                                        exclVal += "product with different melting temperatures detected;"
                            else:
                                if finPeakCount == 1:
                                    if res[oRow][rar_sample_type] in ["ntc", "nac", "ntp", "nrt"]:
                                        exclVal += "product detected in negative control;"

                        if finPeakCount > 1:
                            if res[oRow][rar_sample_type] in ["ntc", "nac", "ntp", "nrt"]:
                                noteVal += "several products with different melting temperatures detected;"
                            if res[oRow][rar_sample_type] in ["std", "pos", "unkn"]:
                                exclVal += "several products with different melting temperatures detected;"

                        # Write back
                        res[oRow][rar_excl] = re.sub(r'^;|;$', '', exclVal)
                        res[oRow][rar_note] = re.sub(r'^;|;$', '', noteVal)

                        resTable[rowPos][0] = res[oRow][rar_id]
                        resTable[rowPos][1] = res[oRow][rar_well]
                        resTable[rowPos][2] = res[oRow][rar_sample]
                        resTable[rowPos][3] = res[oRow][rar_sample_type]
                        resTable[rowPos][4] = res[oRow][rar_tar]
                        resTable[rowPos][5] = res[oRow][rar_excl]
                        resTable[rowPos][6] = res[oRow][rar_note]
                        resTable[rowPos][7] = res[oRow][rar_exp_melt_temp]

                        oCol = truePeakFinPos[oRow]
                        if oCol >= 0:
                            resTable[rowPos][8] = peakResTemp[oRow][oCol]
                            resTable[rowPos][9] = peakResWidth[oRow][oCol]
                            resTable[rowPos][10] = peakResFluor[oRow][oCol]
                            resTable[rowPos][11] = peakResFluor[oRow][oCol] / peakResSumFuor[oRow]
                            resTable[rowPos][12] = peakResDeltaH[oRow][oCol]
                            resTable[rowPos][13] = peakResDeltaH[oRow][oCol] / peakResSumDeltaH[oRow]
                            resTable[rowPos][14] = peakResH[oRow][oCol]
                            resTable[rowPos][15] = peakResH[oRow][oCol] / peakResSumH[oRow]

                            meanResTemp += peakResTemp[oRow][oCol]
                            meanResPeakWidth += peakResWidth[oRow][oCol]
                            meanResFluor += peakResFluor[oRow][oCol]
                            meanResSumFuor += peakResFluor[oRow][oCol] / peakResSumFuor[oRow]
                            meanResDeltaH += peakResDeltaH[oRow][oCol]
                            meanResSumDeltaH += peakResDeltaH[oRow][oCol] / peakResSumDeltaH[oRow]
                            meanResH += peakResH[oRow][oCol]
                            meanResSumH += peakResH[oRow][oCol] / peakResSumH[oRow]
                            countAddedRows += 1
                        else:
                            resTable[rowPos][9] = 0.0
                            resTable[rowPos][10] = 0.0
                            resTable[rowPos][11] = 0.0
                            resTable[rowPos][12] = 0.0
                            resTable[rowPos][13] = 0.0
                            resTable[rowPos][14] = 0.0
                            resTable[rowPos][15] = 0.0

                if countAddedRows > 0:
                    resTable[averageRow][8] = meanResTemp / countAddedRows
                    resTable[averageRow][9] = meanResPeakWidth / countAddedRows
                    resTable[averageRow][10] = meanResFluor / countAddedRows
                    resTable[averageRow][11] = meanResSumFuor / countAddedRows
                    resTable[averageRow][12] = meanResDeltaH / countAddedRows
                    resTable[averageRow][13] = meanResSumDeltaH / countAddedRows
                    resTable[averageRow][14] = meanResH / countAddedRows
                    resTable[averageRow][15] = meanResSumH / countAddedRows
                else:
                    resTable[averageRow][9] = 0.0
                    resTable[averageRow][10] = 0.0
                    resTable[averageRow][11] = 0.0
                    resTable[averageRow][12] = 0.0
                    resTable[averageRow][13] = 0.0
                    resTable[averageRow][14] = 0.0
                    resTable[averageRow][15] = 0.0
                resTable[averageRow][1] = rowPos - averageRow

                # Add all the other peaks
                sortPeaks = sorted(artifactPeaks[curTarNr], key=float)
                workPeaks = []
                lastTemp = 120.0
                truePeakRow = -1
                for curPeak in range(0, len(sortPeaks)):
                    if truePeakRow == -1 and sortPeaks[curPeak] > float(expTemp[curTarNr]) > 0.0:
                        truePeakRow = curPeak
                        workPeaks.append(float(expTemp[curTarNr]))
                    lastTemp = sortPeaks[curPeak]
                    workPeaks.append(lastTemp)
                if truePeakRow == -1 and float(expTemp[curTarNr]) > lastTemp:
                    truePeakRow = len(sortPeaks)
                    workPeaks.append(float(expTemp[curTarNr]))

                for curCol in range(0, len(workPeaks)):
                    rowPos = averageRow
                    if curCol != truePeakRow:
                        artiCol = artifactPeaks[curTarNr].index(workPeaks[curCol])
                    else:
                        artiCol = -1
                    colOffset = 17 + 9 * curCol
                    meanResTemp = 0.0
                    meanResPeakWidth = 0.0
                    meanResH = 0.0
                    meanResSumH = 0.0
                    meanResDeltaH = 0.0
                    meanResSumDeltaH = 0.0
                    meanResFluor = 0.0
                    meanResSumFuor = 0.0
                    countAddedRows = 0
                    for oRow in range(0, spFl[0]):
                        if curTarNr == tarWinLookup[res[oRow][rar_tar]]:
                            rowPos += 1
                            if curCol == truePeakRow:
                                oCol = truePeakFinPos[oRow]
                            else:
                                oCol = orderPeaks[oRow][artiCol]

                            if oCol >= 0:
                                resTable[rowPos][colOffset] = peakResTemp[oRow][oCol]
                                resTable[rowPos][colOffset + 1] = peakResWidth[oRow][oCol]
                                resTable[rowPos][colOffset + 2] = peakResFluor[oRow][oCol]
                                resTable[rowPos][colOffset + 3] = peakResFluor[oRow][oCol] / peakResSumFuor[oRow]
                                resTable[rowPos][colOffset + 4] = peakResDeltaH[oRow][oCol]
                                resTable[rowPos][colOffset + 5] = peakResDeltaH[oRow][oCol] / peakResSumDeltaH[oRow]
                                resTable[rowPos][colOffset + 6] = peakResH[oRow][oCol]
                                resTable[rowPos][colOffset + 7] = peakResH[oRow][oCol] / peakResSumH[oRow]
                                meanResTemp += peakResTemp[oRow][oCol]
                                meanResPeakWidth += peakResWidth[oRow][oCol]
                                meanResFluor += peakResFluor[oRow][oCol]
                                meanResSumFuor += peakResFluor[oRow][oCol] / peakResSumFuor[oRow]
                                meanResDeltaH += peakResDeltaH[oRow][oCol]
                                meanResSumDeltaH += peakResDeltaH[oRow][oCol] / peakResSumDeltaH[oRow]
                                meanResH += peakResH[oRow][oCol]
                                meanResSumH += peakResH[oRow][oCol] / peakResSumH[oRow]
                                countAddedRows += 1
                            else:
                                resTable[rowPos][colOffset + 1] = 0.0
                                resTable[rowPos][colOffset + 2] = 0.0
                                resTable[rowPos][colOffset + 3] = 0.0
                                resTable[rowPos][colOffset + 4] = 0.0
                                resTable[rowPos][colOffset + 5] = 0.0
                                resTable[rowPos][colOffset + 6] = 0.0
                                resTable[rowPos][colOffset + 7] = 0.0
                    if countAddedRows > 0:
                        resTable[averageRow][colOffset - 1] = countAddedRows
                        resTable[averageRow][colOffset] = meanResTemp / countAddedRows
                        resTable[averageRow][colOffset + 1] = meanResPeakWidth / countAddedRows
                        resTable[averageRow][colOffset + 2] = meanResFluor / countAddedRows
                        resTable[averageRow][colOffset + 3] = meanResSumFuor / countAddedRows
                        resTable[averageRow][colOffset + 4] = meanResDeltaH / countAddedRows
                        resTable[averageRow][colOffset + 5] = meanResSumDeltaH / countAddedRows
                        resTable[averageRow][colOffset + 6] = meanResH / countAddedRows
                        resTable[averageRow][colOffset + 7] = meanResSumH / countAddedRows
                    else:
                        resTable[averageRow][colOffset - 1] = 0
                        resTable[averageRow][colOffset + 1] = 0.0
                        resTable[averageRow][colOffset + 2] = 0.0
                        resTable[averageRow][colOffset + 3] = 0.0
                        resTable[averageRow][colOffset + 4] = 0.0
                        resTable[averageRow][colOffset + 5] = 0.0
                        resTable[averageRow][colOffset + 6] = 0.0
                        resTable[averageRow][colOffset + 7] = 0.0

            ##############################
            # write out the rdml results #
            ##############################
            if updateRDML is True:
                dataXMLelements = _getXMLDataType()
                for rRow in range(0, len(res)):
                    if rdmlElemData[rRow] is not None:
                        _change_subelement(rdmlElemData[rRow], "excl", dataXMLelements, res[rRow][rar_excl], True, "string")
                        if dataVersion in ["1.3", "1.4"]:
                            _change_subelement(rdmlElemData[rRow], "note", dataXMLelements, res[rRow][rar_note], True, "string")
                            lCol = truePeakFinPos[rRow]
                            if lCol >= 0:
                                if res[rRow][rar_tar_chemistry] == "saturating DNA binding dye":
                                    finalFactor = peakResFluor[rRow][lCol] / peakResSumFuor[rRow]
                                    goodVal = "{:.3f}".format(finalFactor)
                                    _change_subelement(rdmlElemData[rRow], "corrF", dataXMLelements, goodVal, True, "string")
                                    _change_subelement(rdmlElemData[rRow], "corrCq", dataXMLelements, "", True, "string")
                            else:
                                if res[rRow][rar_tar_chemistry] == "saturating DNA binding dye":
                                     _change_subelement(rdmlElemData[rRow], "corrF", dataXMLelements, "0.0", True, "string")
                                     _change_subelement(rdmlElemData[rRow], "corrCq", dataXMLelements, "", True, "string")
            finalData["resultsList"] = resTable

        if err != "":
            finalData["error"] = err

        return finalData


def main():
    parser = argparse.ArgumentParser(description='The command line interface to the RDML-Python library.')
    parser.add_argument('--version', action='store_true', help='print version number')
    parser.add_argument('-v', '--validate', metavar="data.rdml", dest='validate', help='validate file against schema')
    parser.add_argument('-e', '--experiment', metavar="exp_1", dest='experiment', help='select experiment')
    parser.add_argument('-le', '--listexperiment', metavar="data.rdml", dest='listExp', help='list experiments')
    parser.add_argument('-r', '--run', metavar="run_1", dest='run', help='select run')
    parser.add_argument('-lr', '--listrun', metavar="data.rdml", dest='listRun', help='list runs')
    parser.add_argument('-lrp', '--linRegPCR', metavar="data.rdml", dest='linRegPCR', help='run LinRegPCR')
    parser.add_argument('-o', '--resultfile', metavar="data_out.rdml", dest='resultfile',
                        help='LinRegPCR: output file of LinRegPCR')
    parser.add_argument('--pcrEfficiencyExl', metavar="0.05",
                        help='LinRegPCR: provide a range for for exclusion from mean PCR efficiency')
    parser.add_argument('--excludeNoPlateau', action='store_true',
                        help='LinRegPCR: exclude no plateau samples from mean PCR efficiency')
    parser.add_argument('--includeInstableBaseline', action='store_true',
                        help='LinRegPCR: include samples with an instable baseline')
    parser.add_argument('--excludeEfficiency', metavar="outlier",
                        help='LinRegPCR: choose [outlier, mean, include] to exclude different individual efficiency ' +
                             'samples from mean PCR efficiency')
    parser.add_argument('--commaConv', action='store_true', help='LinRegPCR: convert comma to dot in numbers')
    parser.add_argument('--ignoreExclusion', action='store_true', help='LinRegPCR: ignore the exclusion field')
    parser.add_argument('--saveRaw', metavar="raw_data.csv",
                        help='LinRegPCR & MeltCurveAnalysis: output file for raw (unmodified) data')
    parser.add_argument('--saveBaslineCorr', metavar="baseline_corrected.csv",
                        help='LinRegPCR: output file for baseline corrected data')
    parser.add_argument('-mca', '--meltCurveAnalysis', metavar="data.rdml", dest='meltCurveAnalysis',
                        help='run MeltCurveAnalysis')
    parser.add_argument('--mcaNormMethod', metavar="exponential",
                        help='MeltCurveAnalysis: choose [exponential, bilinear, combined] for baseline correction')
    parser.add_argument('--mcaFluorSource', metavar="normalised",
                        help='MeltCurveAnalysis: choose [normalised, smooth] for fluorescence correction factor calculation')
    parser.add_argument('--mcaTruePeakWidth', metavar="1.0",
                        help='MeltCurveAnalysis: provide the maximum allowed width in Celsius of the expected (true) peak')
    parser.add_argument('--mcaArtifactPeakWidth', metavar="1.0",
                        help='MeltCurveAnalysis: provide the maximum allowed width in Celsius of all artifact peaks')
    parser.add_argument('--mcaExpoLowTemp', metavar="65.0",
                        help='MeltCurveAnalysis: provide the low temperature for the exponential normalisation')
    parser.add_argument('--mcaExpoHighTemp', metavar="92.0",
                        help='MeltCurveAnalysis: provide the high temperature for the exponential normalisation')
    parser.add_argument('--mcaBilinLowStartTemp', metavar="68.0",
                        help='MeltCurveAnalysis: provide the low start temperature for the bilinear normalisation')
    parser.add_argument('--mcaBilinLowStopTemp', metavar="70.0",
                        help='MeltCurveAnalysis: provide the low stop temperature for the bilinear normalisation')
    parser.add_argument('--mcaBilinHighStartTemp', metavar="93.0",
                        help='MeltCurveAnalysis: provide the high start temperature for the bilinear normalisation')
    parser.add_argument('--mcaBilinHighStopTemp', metavar="94.0",
                        help='MeltCurveAnalysis: provide the high stop temperature for the bilinear normalisation')
    parser.add_argument('--mcaPeakLowTemp', metavar="60.0",
                        help='mcaPeakLowTemp: peaks below this temperature are ignored')
    parser.add_argument('--mcaPeakHighTemp', metavar="98.0",
                        help='mcaPeakHighTemp: peaks above this temperature are ignored')
    parser.add_argument('--mcaPeakMaxWidth', metavar="5.0",
                        help='mcaPeakMaxWidth: peaks broader than this temperature width are ignored')
    parser.add_argument('--mcaPeakCutoff', metavar="5.0",
                        help='mcaPeakCutoff: the percentage below melting peaks are ignored in calculations')
    parser.add_argument('--saveDerivative', metavar="processed_data",
                        help='MeltCurveAnalysis: base name for calculated derivative data')
    parser.add_argument('--saveResults', metavar="results.csv", help='LinRegPCR & MeltCurveAnalysis: output results as table')
    parser.add_argument('--timeRun', action='store_true', help='LinRegPCR: print a timestamp')
    parser.add_argument('--verbose', action='store_true', help='LinRegPCR & MeltCurveAnalysis: print comments')

    parser.add_argument("-d", "--doooo", dest="doooo", help="just do stuff")

    args = parser.parse_args()

    if args.version:
        print("rdmlpython version " + get_rdml_lib_version())
        sys.exit(0)

    # Validate RDML file
    if args.validate:
        cli_validate = Rdml()
        cli_resValidate = cli_validate.validate(filename=args.validate)
        print(cli_resValidate)
        sys.exit(0)

    # List all Experiments
    if args.listExp:
        cli_listExp = Rdml(args.listExp)
        cli_expList = cli_listExp.experiments()
        print("Experiments in file \"" + args.listExp + "\":")
        if len(cli_expList) < 1:
            print("No experiments found!")
            sys.exit(0)
        for cli_exp in cli_expList:
            print(cli_exp["id"])
        sys.exit(0)

    # List all Runs
    if args.listRun:
        cli_listRun = Rdml(args.listRun)
        if args.experiment:
            try:
                cli_exp = cli_listRun.get_experiment(byid=args.experiment)
            except RdmlError as cli_err:
                print("Error: " + str(cli_err))
                sys.exit(1)
            else:
                print("Using experiment: \"" + args.experiment + "\"")
        else:
            cli_expList = cli_listRun.experiments()
            if len(cli_expList) < 1:
                print("No experiments found!")
                sys.exit(0)
            cli_exp = cli_expList[0]
            print("No experiment given (use option -e). Using \"" + cli_expList[0]["id"] + "\"")

        cli_runList = cli_exp.runs()
        print("Runs in file \"" + args.listRun + "\":")
        if len(cli_runList) < 1:
            print("No runs found!")
            sys.exit(0)
        for cli_run in cli_runList:
            print(cli_run["id"])
        sys.exit(0)

    # Run LinRegPCR from commandline
    if args.linRegPCR:
        cli_linRegPCR = Rdml(args.linRegPCR)
        if cli_linRegPCR.version() == "1.0":
            cli_linRegPCR.migrate_version_1_0_to_1_1()
        if args.experiment:
            try:
                cli_exp = cli_linRegPCR.get_experiment(byid=args.experiment)
            except RdmlError as cli_err:
                print("Error: " + str(cli_err))
                sys.exit(1)
            else:
                print("Using experiment: \"" + args.experiment + "\"")
        else:
            cli_expList = cli_linRegPCR.experiments()
            if len(cli_expList) < 1:
                print("No experiments found!")
                sys.exit(0)
            cli_exp = cli_expList[0]
            print("No experiment given (use option -e). Using \"" + cli_expList[0]["id"] + "\"")
        if args.run:
            try:
                cli_run = cli_exp.get_run(byid=args.run)
            except RdmlError as cli_err:
                print("Error: " + str(cli_err))
                sys.exit(1)
            else:
                print("Using run: \"" + args.run + "\"")
        else:
            cli_runList = cli_exp.runs()
            if len(cli_runList) < 1:
                print("No runs found!")
                sys.exit(0)
            cli_run = cli_runList[0]
            print("No run given (use option -r). Using \"" + cli_runList[0]["id"] + "\"")

        cli_pcrEfficiencyExl = 0.05
        cli_excludeNoPlateau = False
        cli_excludeEfficiency = "outlier"
        cli_excludeInstableBaseline = True
        cli_commaConv = False
        cli_ignoreExclusion = False
        cli_timeRun = False
        cli_verbose = False
        cli_saveRDML = False
        cli_saveRawData = False
        cli_saveBaselineData = False
        cli_saveResultData = False

        if args.pcrEfficiencyExl:
            cli_pcrEfficiencyExl = float(args.pcrEfficiencyExl)
        if args.excludeNoPlateau:
            cli_excludeNoPlateau = True
        if args.excludeEfficiency:
            cli_excludeEfficiency = args.excludeEfficiency
        if args.includeInstableBaseline:
            cli_excludeInstableBaseline = False
        else:
            cli_excludeInstableBaseline = True
        if args.commaConv:
            cli_commaConv = True
        if args.ignoreExclusion:
            cli_ignoreExclusion = True
        if args.timeRun:
            cli_timeRun = True
        if args.verbose:
            cli_verbose = True
        if args.resultfile:
            cli_saveRDML = True
        if args.saveRaw:
            cli_saveRawData = True
        if args.saveBaslineCorr:
            cli_saveBaselineData = True
        if args.saveResults:
            cli_saveResultData = True

        cli_result = cli_run.linRegPCR(pcrEfficiencyExl=cli_pcrEfficiencyExl, updateRDML=cli_saveRDML,
                                       excludeNoPlateau=cli_excludeNoPlateau, excludeEfficiency=cli_excludeEfficiency,
                                       excludeInstableBaseline=cli_excludeInstableBaseline,
                                       commaConv=cli_commaConv, ignoreExclusion=cli_ignoreExclusion,
                                       saveRaw=cli_saveRawData, saveBaslineCorr=cli_saveBaselineData,
                                       saveResultsList=False, saveResultsCSV=cli_saveResultData,
                                       timeRun=cli_timeRun, verbose=cli_verbose)

        if "error" in cli_result:
            print(cli_result["error"])
        if args.resultfile:
            cli_linRegPCR.save(args.resultfile)
        if args.saveRaw:
            if "rawData" in cli_result:
                with open(args.saveRaw, "w") as cli_f:
                    cli_ResStr = ""
                    for cli_row in cli_result["rawData"]:
                        for cli_col in cli_row:
                            if type(cli_col) is float:
                                cli_ResStr += "{0:0.3f}".format(cli_col) + "\t"
                            else:
                                cli_ResStr += str(cli_col) + "\t"
                        cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                    cli_f.write(cli_ResStr)
        if args.saveBaslineCorr:
            if "baselineCorrectedData" in cli_result:
                with open(args.saveBaslineCorr, "w") as cli_f:
                    cli_ResStr = ""
                    for cli_row in cli_result["baselineCorrectedData"]:
                        for cli_col in cli_row:
                            if type(cli_col) is float:
                                cli_ResStr += "{0:0.6f}".format(cli_col) + "\t"
                            else:
                                cli_ResStr += str(cli_col) + "\t"
                        cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                    cli_f.write(cli_ResStr)
        if args.saveResults:
            with open(args.saveResults, "w") as cli_f:
                cli_f.write(cli_result["resultsCSV"])
        sys.exit(0)

    # Run meltCurveAnalysis from commandline
    if args.meltCurveAnalysis:
        cli_meltCurveAnalysis = Rdml(args.meltCurveAnalysis)
        if cli_meltCurveAnalysis.version() == "1.0":
            cli_meltCurveAnalysis.migrate_version_1_0_to_1_1()
        if args.experiment:
            try:
                cli_exp = cli_meltCurveAnalysis.get_experiment(byid=args.experiment)
            except RdmlError as cli_err:
                print("Error: " + str(cli_err))
                sys.exit(1)
            else:
                print("Using experiment: \"" + args.experiment + "\"")
        else:
            cli_expList = cli_meltCurveAnalysis.experiments()
            if len(cli_expList) < 1:
                print("No experiments found!")
                sys.exit(0)
            cli_exp = cli_expList[0]
            print("No experiment given (use option -e). Using \"" + cli_expList[0]["id"] + "\"")
        if args.run:
            try:
                cli_run = cli_exp.get_run(byid=args.run)
            except RdmlError as cli_err:
                print("Error: " + str(cli_err))
                sys.exit(1)
            else:
                print("Using run: \"" + args.run + "\"")
        else:
            cli_runList = cli_exp.runs()
            if len(cli_runList) < 1:
                print("No runs found!")
                sys.exit(0)
            cli_run = cli_runList[0]
            print("No run given (use option -r). Using \"" + cli_runList[0]["id"] + "\"")

        cli_saveRDML = False
        cli_saveRawData = False
        cli_saveDerivative = False
        cli_saveResultData = False
        cli_mcaNormMethod = "exponential"
        cli_mcaFluorSource = "normalised"
        cli_mcaTruePeakWidth = 1.0
        cli_mcaArtifactPeakWidth = 1.0
        cli_mcaExpoLowTemp = 65.0
        cli_mcaExpoHighTemp = 92.0
        cli_mcaBilinLowStartTemp = 65.0
        cli_mcaBilinLowStopTemp = 67.0
        cli_mcaBilinHighStartTemp = 93.0
        cli_mcaBilinHighStopTemp = 94.0
        cli_mcaPeakLowTemp = 60.0
        cli_mcaPeakHighTemp = 98.0
        cli_mcaPeakMaxWidth = 5.0
        cli_mcaPeakCutoff = 5.0

        if args.resultfile:
            cli_saveRDML = True
        if args.saveRaw:
            cli_saveRawData = True
        if args.saveDerivative:
            cli_saveDerivative = True
        if args.saveResults:
            cli_saveResultData = True
        if args.mcaNormMethod:
            cli_mcaNormMethod = args.mcaNormMethod
        if args.mcaFluorSource:
            cli_mcaFluorSource = args.mcaFluorSource
        if args.mcaTruePeakWidth:
            cli_mcaTruePeakWidth = float(args.mcaTruePeakWidth)
        if args.mcaArtifactPeakWidth:
            cli_mcaArtifactPeakWidth = float(args.mcaArtifactPeakWidth)
        if args.mcaExpoLowTemp:
            cli_mcaExpoLowTemp = float(args.mcaExpoLowTemp)
        if args.mcaExpoHighTemp:
            cli_mcaExpoHighTemp = float(args.mcaExpoHighTemp)
        if args.mcaBilinLowStartTemp:
            cli_mcaBilinLowStartTemp = float(args.mcaBilinLowStartTemp)
        if args.mcaBilinLowStopTemp:
            cli_mcaBilinLowStopTemp = float(args.mcaBilinLowStopTemp)
        if args.mcaBilinHighStartTemp:
            cli_mcaBilinHighStartTemp = float(args.mcaBilinHighStartTemp)
        if args.mcaBilinHighStopTemp:
            cli_mcaBilinHighStopTemp = float(args.mcaBilinHighStopTemp)
        if args.mcaPeakLowTemp:
            cli_mcaPeakLowTemp = float(args.mcaPeakLowTemp)
        if args.mcaPeakHighTemp:
            cli_mcaPeakHighTemp = float(args.mcaPeakHighTemp)
        if args.mcaPeakMaxWidth:
            cli_mcaPeakMaxWidth = float(args.mcaPeakMaxWidth)
        if args.mcaPeakCutoff:
            cli_mcaPeakCutoff = float(args.mcaPeakCutoff)

        cli_result = cli_run.meltCurveAnalysis(normMethod=cli_mcaNormMethod,
                                               fluorSource=cli_mcaFluorSource,
                                               truePeakWidth=cli_mcaTruePeakWidth,
                                               artifactPeakWidth=cli_mcaArtifactPeakWidth,
                                               expoLowTemp=cli_mcaExpoLowTemp,
                                               expoHighTemp=cli_mcaExpoHighTemp,
                                               bilinLowStartTemp=cli_mcaBilinLowStartTemp,
                                               bilinLowStopTemp=cli_mcaBilinLowStopTemp,
                                               bilinHighStartTemp=cli_mcaBilinHighStartTemp,
                                               bilinHighStopTemp=cli_mcaBilinHighStopTemp,
                                               peakLowTemp=cli_mcaPeakLowTemp,
                                               peakHighTemp=cli_mcaPeakHighTemp,
                                               peakMaxWidth=cli_mcaPeakMaxWidth,
                                               peakCutoff=cli_mcaPeakCutoff,
                                               updateRDML=cli_saveRDML,
                                               saveRaw=cli_saveRawData,
                                               saveDerivative=cli_saveDerivative,
                                               saveResultsList=True,
                                               saveResultsCSV=cli_saveResultData)
        if args.saveRaw:
            if "rawData" in cli_result:
                with open(args.saveRaw, "w") as cli_f:
                    cli_ResStr = ""
                    for cli_row in cli_result["rawData"]:
                        for cli_col in cli_row:
                            if type(cli_col) is float:
                                cli_ResStr += "{0:0.3f}".format(float(cli_col)) + "\t"
                            else:
                                cli_ResStr += str(cli_col) + "\t"
                        cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                    cli_f.write(cli_ResStr)

        if args.saveDerivative:
            if "derivative" in cli_result:
                if "smoothed" in cli_result["derivative"]:
                    with open(args.saveDerivative + "_smoothed.tsv", "w") as cli_f:
                        cli_ResStr = ""
                        for cli_row in cli_result["derivative"]["smoothed"]:
                            for cli_col in cli_row:
                                if type(cli_col) is float:
                                    cli_ResStr += "{0:0.8f}".format(float(cli_col)) + "\t"
                                else:
                                    cli_ResStr += str(cli_col) + "\t"
                            cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                        cli_f.write(cli_ResStr)

                if "normalized" in cli_result["derivative"]:
                    with open(args.saveDerivative + "_normalized.tsv", "w") as cli_f:
                        cli_ResStr = ""
                        for cli_row in cli_result["derivative"]["normalized"]:
                            for cli_col in cli_row:
                                if type(cli_col) is float:
                                    cli_ResStr += "{0:0.8f}".format(float(cli_col)) + "\t"
                                else:
                                    cli_ResStr += str(cli_col) + "\t"
                            cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                        cli_f.write(cli_ResStr)

                if "firstDerivative" in cli_result["derivative"]:
                    with open(args.saveDerivative + "_firstDerivative.tsv", "w") as cli_f:
                        cli_ResStr = ""
                        for cli_row in cli_result["derivative"]["firstDerivative"]:
                            for cli_col in cli_row:
                                if type(cli_col) is float:
                                    cli_ResStr += "{0:0.8e}".format(float(cli_col)) + "\t"
                                else:
                                    cli_ResStr += str(cli_col) + "\t"
                            cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                        cli_f.write(cli_ResStr)

                if "secondDerivative" in cli_result["derivative"]:
                    with open(args.saveDerivative + "_secondDerivative.tsv", "w") as cli_f:
                        cli_ResStr = ""
                        for cli_row in cli_result["derivative"]["secondDerivative"]:
                            for cli_col in cli_row:
                                if type(cli_col) is float:
                                    cli_ResStr += "{0:0.8e}".format(float(cli_col)) + "\t"
                                else:
                                    cli_ResStr += str(cli_col) + "\t"
                            cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                        cli_f.write(cli_ResStr)

        if args.saveResults:
            if "resultsList" in cli_result:
                with open(args.saveResults, "w") as cli_f:
                    cli_ResStr = ""
                    for cli_row in cli_result["resultsList"]:
                        for cli_col in cli_row:
                            if isinstance(cli_col, np.float64) or isinstance(cli_col, float):
                                cli_ResStr += "{0:0.3f}".format(float(cli_col)) + "\t"
                            else:
                                cli_ResStr += str(cli_col) + "\t"
                        cli_ResStr = re.sub(r"\t$", "\n", cli_ResStr)
                    cli_f.write(cli_ResStr)

        sys.exit(0)


if __name__ == "__main__":
    main()
