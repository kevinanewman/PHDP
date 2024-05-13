"""

**Various functions that may be used throughout PHDP.**


----

**CODE**

"""

import sys
import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np

# np.seterr(all='raise')  # for troubleshooting runtime warnings


def get_ip_address():
    """
    Attempt to get "local" IP address(es)

    Example:

    ::

        >>> socket.gethostbyname_ex(socket.gethostname())
        ('mac-mini.local', [], ['127.0.0.1', '192.168.1.20'])

    Returns: list of local IP address(es)

    """
    import socket

    my_ip = []

    retries = 0
    ip_found = False
    while not ip_found and retries < 10:
        try:
            my_ip = socket.gethostbyname_ex(socket.gethostname())[2]
            ip_found = True
        except:
            retries += 1

    if not my_ip.count('127.0.0.1'):
        my_ip.append('127.0.0.1')  # Add support for local loopback interface

    return my_ip


def dataframe_to_numeric(df):
    """
    Convert dataframe columns to numeric (i.e. non-object dtypes) if possible.

    Args:
        df (DataFrame): the dataframe to convert to numeric

    Returns:
        df with numeric columns where possible

    """
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='ignore')

    return df


def series_to_numeric(ser):
    """
    Convert series entries to numeric (i.e. non-object dtypes) if possible.

    Args:
        ser (Series): the series to convert to numeric

    Returns:
        ser with numeric columns where possible

    """
    ser_out = pd.Series(dtype='float64')
    for c in ser.keys():
        ser_out[c] = pd.to_numeric(ser[c], errors='ignore')

    return ser_out


def sum_dict(dict_in, include=None, exclude=None):
    """
    Add up all terms in a dict given the ``include`` and ``exclude`` constraints.

    Args:
        dict_in (numeric dict_like): the object with elements to sum
        include (str): include term in sum if ``include`` in dict key
        exclude (str): exclude term from some if ``exclude`` in dict key

    Returns:
        Sum of terms in ``dict_in`` given the include and exclude constraints.

    """
    keys = sorted(dict_in.keys())
    if include is not None:
        keys = [k for k in keys if include in k]
    if exclude is not None:
        keys = [k for k in keys if exclude not in k]

    return sum([dict_in[k] for k in keys])


def print_keys(dict_in, include=None, exclude=None, values=True):
    """
    Print some or all keys (and optionally values) in a dict-like object

    Args:
        dict_in (dict-like): the object with keys to print
        include (str): a substring that must be present, if provided
        exclude (str): a substring that must not be present, if provided
        values (bool): print values if ``True``

    """
    keys = sorted(dict_in.keys())
    if include is not None:
        keys = [k for k in keys if include in k]
    if exclude is not None:
        keys = [k for k in keys if exclude not in k]

    max_key_len = max([len(k) for k in keys])

    for k in keys:
        if values:
            format_str = '%' + '%d' % max_key_len + 's: %s'
            print(format_str % (k, dict_in[k]))
        else:
            print(k)

    return keys


def print_dict(dict_in, num_tabs=0, to_string=False):
    """
    Attempt to printy-print a dictionary to the Python console.

    Args:
        dict_in (dict): dictionary to print
        num_tabs (int): optional argument, used to indent subsequent layers of the dictionary
        to_string (Bool): if True then result will be returned as a printable string, instead of printed to the console

    Returns:
        print_dict string if to_string==True

    """
    s = ''

    if num_tabs == 0:
        s += '\n'
        if type(dict_in) is pd.Series:
            dict_in = dict_in.to_dict()

    if type(dict_in) is list or type(dict_in) is not dict:
        s += '\t' * num_tabs + str(dict_in) + '\n'
    else:
        try:
            keys = sorted(dict_in.keys())
        except:
            keys = dict_in.keys()

        for k in keys:
            if type(dict_in[k]) == list:
                if dict_in[k]:
                    s += '\t' * num_tabs + str(k) + ':' + str(dict_in[k]) + '\n'
                else:
                    s += '\t' * num_tabs + str(k) + '\n'
            else:
                s += '\t' * num_tabs + str(k) + '\n'
                s += print_dict(dict_in[k], num_tabs + 1, to_string=True)

    if num_tabs == 0:
        s += '\n'

    if not to_string:
        print(s[:-1])

    if to_string:
        return s


def print_list(list_in):
    """
    Print the given list, one line per item

    Args:
        list_in (list): the list to print

    """
    for i in list_in:
        print(i)
    print()


def linspace(min_val, max_val, num_values):
    """
    Create a list of num_values evenly spaced values between min and max.  Based on ``Matlab`` linspace command.

    Args:
        min_val (numeric): the minimum value
        max_val (numeric): the maximum value
        num_values (int): the total number of values to return

    Returns:
        A list of evenly spaced values between min and max

    """
    ans = np.arange(min_val, max_val + (max_val - min_val) / (num_values - 1), (max_val - min_val) / (num_values - 1))
    return ans[0:num_values]


def unique(vector):
    """
    Return unique values in a list of values, in order of appearance.

    Args:
        vector ([numeric]): list of values

    Returns:
        List of unique values, in order of appearance

    """
    indexes = np.unique(vector, return_index=True)[1]
    return [vector[index] for index in sorted(indexes)]


def ASTM_round(var, precision=0):
    """
    Rounds numbers as defined in ISO / IEC / IEEE 60559

    Args:
        var (float, Series): number to be rounded, scalar or pandas Series
        precision (int): number of decimal places in result

    Returns:
        var rounded using ASTM method with precision decimal places in result

    """

    if type(var) is list:
        var = np.array(var)

    scaled_var = var * (10 ** precision)

    z = np.remainder(scaled_var, 2)

    if type(z) is pd.core.series.Series or type(z) is np.ndarray:
        if type(z) is np.ndarray:
            z = pd.Series(z)
        z.loc[z != 0.5] = 0
    else:
        if abs(z) != 0.5:
            z = 0

    rounded_number = np.round(scaled_var - z) / (10**precision)

    return rounded_number


def send_text(dest, message, email, password):
    """

    SMS Gateways for each Carrier
    AT&T: [number]@txt.att.net
    Sprint: [number]@messaging.sprintpcs.com or [number]@pm.sprint.com
    T-Mobile: [number]@tmomail.net
    Verizon: [number]@vtext.com
    Boost Mobile: [number]@myboostmobile.com
    Cricket: [number]@sms.mycricket.com
    Metro PCS: [number]@mymetropcs.com
    Tracfone: [number]@mmst5.tracfone.com
    U.S. Cellular: [number]@email.uscc.net
    Virgin Mobile: [number]@vmobl.com

    Args:
        dest (str): e.g. '8005552323@myboostmobile.com'
        message (str): the message to send
        email (str): the email address of the email server to use, e.g. 'my_email@gmail.com'
        password (str): the password for the email account, recommend setting up an app-specific password

    """
    import time
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    sms_gateway = dest

    pas = password

    # The server we use to send emails in our case it will be gmail but every email provider has a different smtp
    # and port is also provided by the email provider.
    smtp = "smtp.gmail.com"
    port = 587
    # This will start our email server
    server = smtplib.SMTP(smtp, port)
    # Starting the server
    server.starttls()
    # Now we need to login
    server.login(email, pas)

    # Now we use the MIME module to structure our message.
    msg = MIMEMultipart()
    msg['From'] = email
    msg['To'] = sms_gateway

    # Make sure you add a new line in the subject
    timestamp_str = time.strftime('%m/%d/%Y %H:%M:%S')
    msg['Subject'] = "%s" % timestamp_str + "\n"

    # Make sure you also add new lines to your body
    if not message.endswith('\n'):
        message = message + '\n'

    # and then attach that body furthermore you can also send html content.
    msg.attach(MIMEText(message, 'plain'))
    sms = msg.as_string()
    server.sendmail(email, sms_gateway, sms)

    # lastly quit the server
    server.quit()


import string
num2alphadict = dict(zip(range(0, 26), string.ascii_uppercase))


def get_letters_from_index(index):
    """
    Produces a string from numbers, like Excel column labels

    0 -> A
    1 -> B
    ...
    26 -> AA
    27 -> AB
    etc ...

    Args:
        index (int): e.g. column number

    Returns:
        Alphabetic string sequence based on num

    """

    alpha_str = ""
    numloops = (index) // 26

    if numloops > 0:
        alpha_str = alpha_str + get_letters_from_index(numloops-1)

    remainder = index % 26

    if remainder > 0:
        alpha_str = alpha_str + num2alphadict[remainder]
    else:
        alpha_str = alpha_str + "A"

    return alpha_str


def is_container(obj):
    """
    Checks if the given object is a non-string iterable, e.g. a container.

    Parameters:
        obj : Any type
            The object to be checked.

    Returns:
        bool: True if the object is iterable, False otherwise.
    """
    from collections.abc import Iterable
    return isinstance(obj, Iterable) and not isinstance(obj, str)


def find_cell(df: pd.DataFrame, search_str: str) -> tuple:
    """
    Finds the row and column number of a cell in a DataFrame that equals a given string.
    Parameters:
        df (pd.DataFrame): The DataFrame to search within.
        search_str (str): The string to search for in the DataFrame.
    Returns:
        tuple: A tuple containing the row and column number of the cell, or ``(None, None)`` if the string is not found.
               If the string is found multiple times, returns the row and column numbers of the first occurrence.
    """
    search_result = (df == search_str).values.nonzero()

    if len(search_result[0]):
        row, col = search_result
        return row[0], col[0]  # return location of first instance
    else:
        return None, None


def set_value_at(df, search_str, set_value, col_offset=1, small_float_format='%.4f', large_float_format='%.3f',
                 large_cutoff=100):
    """
    Finds a given string in a DataFrame and set the value of the cell(s) next to it, as determined by ``col_offset``

    Parameters:
        df (DataFrame): The input DataFrame.
        search_str (str): The string to be searched for
        set_value (obj, iterable): the value(s) to set next to the search_str, if found
        col_offset (int): the column offset from the search target cell to the target value cell
        small_float_format (str): format string for "small" floating point numbers
        large_float_format (str): format string for "large" floating point numbers
        large_cutoff (numeric): numeric cutoff for using large float format

    Returns:
        Raises an Exception if the string is not found, otherwise sets appropriate cell value(s) in ``df``

    """

    row_index, col_index = find_cell(df, search_str)

    if row_index is None:
        raise Exception('"%s" not found in dataframe' % search_str)
    else:
        if is_container(set_value):
            for idx, v in enumerate(set_value):
                if v is not None and type(v) is not str:
                    if v >= large_cutoff:
                        df.iloc[row_index, col_index + col_offset + idx] = large_float_format % v
                    else:
                        df.iloc[row_index, col_index + col_offset + idx] = small_float_format % v
                else:
                    df.iloc[row_index, col_index + col_offset + idx] = v
        else:
            df.iloc[row_index, col_index + col_offset] = set_value


def pass_fail_range(value, allowed_range):
    """
    This function determines whether a given value falls within a specified range.
    It takes two arguments: the value to be checked (`value`) and the range of acceptable values (`range`).

    Arguments:
        value (any): The value to be evaluated in relation to the provided range.
        allowed_range (list or tuple of 2 elements): A list or tuple consisting of the lower and upper bounds of the valid range for the value.

    Returns:
        str: Returns 'pass' if the `value` falls within the range, otherwise returns 'FAIL'.

    """
    if allowed_range[0] <= value <= allowed_range[1]:
        return 'pass'
    else:
        return 'FAIL'


def handle_emscal_driftline(signal_name):
    """
    Get driftline and component name from signal, to look up values in EmsCaResults

    Args:
        signal_name (str): signal name, suchas 'conRawCO2_Avg_%vol'

    Returns:
        EmsCal component name, driftline and scale_factor

    """
    rename = {
        'LCO': 'COL',
        'HCO': 'COH',
        'CH4cutter': 'CH4'
    }

    signal, signal_type, unit = signal_name.split('_')

    if signal.startswith('conRaw'):
        component = signal.replace('conRaw', '')
        driftline = 'DIRECT'
        if component == 'NH3':
            driftline = 'HOT'
    elif signal.startswith('conEGRCO2'):
        component = 'CO2'
        driftline = 'EGR'
    else:
        component = signal.replace('con', '')
        driftline = 'DILUTE'

    if component in rename:
        component = rename[component]

    if unit == '%vol':
        scale_factor = 10 ** 4
    else:
        scale_factor = 1

    return component, driftline, scale_factor


if __name__ == '__main__':
    try:
        pass
    except:
        import os
        import traceback
        print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        os._exit(-1)
