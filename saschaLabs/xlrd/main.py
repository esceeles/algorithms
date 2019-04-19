# import the packages `xlrd` and `statistics`
import xlrd as x
import statistics as s

def open_fish_file():
    """
    Opens `pdx_fish.xls` and returns contents as column-major list of
    lists. Dates are left as floating point values in Excel's epoch format.
    """
    pdx_fish_file = "./pdx_fish.xls"

    # open workbook, creating an object of the xlrd `Book` class
    # recommended xlrd function `open_workbook`
    book = x.open_workbook(pdx_fish_file)

    # assign the first worksheet in the workbook you just created to the
    # name `pdx_fish` recommended xlrd `Book` method: `sheet_by_index()`
    pdx_fish = book.sheet_by_index(0)

    # the number of columns in the worksheet
    # recommended xlrd `Book` attribute : `ncols`
    quant_columns = pdx_fish.ncols

    # the number of rows in the worksheet
    # recommended xlrd `Book` attribute: `nrows`
    quant_rows = pdx_fish.nrows

    # iterates through xlrd worksheet and copies values to list of lists
    worksheet = list()
    for col in range(quant_columns):
        worksheet.append(list())
        for row in range(quant_rows):
            worksheet[col].append(pdx_fish.cell_value(row, col))

    return worksheet


def print_worksheet_xlrd(worksheet_xlrd):
    """
    Takes a worksheet in the form of an xlrd worksheet object
    and prints each column to std out.
    """
    for col in range(len(worksheet_xlrd)):
        for row in range(len(worksheet_xlrd[0])):
            print(worksheet_xlrd[col][row])
        print("\n")


def print_worksheet_lol(worksheet_lol):
    """
    Takes a worksheet in the form of a column-major list of lists
    and prints each column to std out.
    """
    for col in range(len(worksheet_lol)):
        for row in range(len(worksheet_lol[0])):
            print(worksheet_lol[col][row])
        print("\n")


def header_col_number_mapping(worksheet_lol):
    """
    Takes a worksheet in the form of a column-major list of lists
    and returns a dictionary mapping the names of the columns to the
    column numbers, indexed from zero.
    """
    header_col_number_dict = dict()
    for col in range(len(worksheet_lol)):
        header_col_number_dict[worksheet_lol[col][0]] = col
    return header_col_number_dict


def nums_only_in_list(lst):
    """
    Takes a list of mixed data types and returns a list containing
    only the `int` and `float` values in the original list.
    """
    return [num for num in lst if isinstance(num, (int, float))]


if __name__ == "__main__":
    wb = open_fish_file()
    names = header_col_number_mapping(wb)
    # Using the `statistics` package, find the mean and standard deviation
    # of the lengths of the fish in `pdx_fish.xls`, disregarding all
    # non-numeric values.
    col12 = wb[12]
    lengths = nums_only_in_list(col12)
    print("mean: {}".format(s.mean(lengths)))
    print("standard deviation: " + str(s.stdev(lengths)))
    # From the file `pdx_fish.xls`, what is the largest Reticulate Sculpin
    # obsverved? (again, disregarding non-numeric values)
    rs = list()
    count = 0
    for i in wb[10]:
      if "Reticulate Sculpin" in i:
         rs.append(wb[12][count])
      count += 1
    rsnums = nums_only_in_list(rs)
    print("largest Reticulate Sculpin: " + str(max(rsnums)))
