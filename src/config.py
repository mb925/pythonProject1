import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/../'

data = {
    "backup": absolute + 'data/backup',
    "mobi": absolute + 'data/mobi',
    "repeatsdb": absolute + 'data/repeatsdblite/repeatsdb'
}
