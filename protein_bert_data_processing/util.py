import sys
import os
from datetime import datetime, timedelta

### Logging ###

def log(*message, **kwargs):
    
    global _log_file
    
    end = kwargs.get('end', '\n')
    
    if len(message) == 1:
        message, = message
    
    full_message = '[%s] %s' % (format_now(), message)
    
    print(full_message, end = end)
    sys.stdout.flush()
    
    if log_file_open():
        _log_file.write(full_message + end)
        _log_file.flush()


def start_log(log_dir, log_file_base_name):
    
    global _log_file
    
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    log_file_name = '%s__%d__%s.txt' % (log_file_base_name, os.getpid(), format_now())
    
    if not log_file_open():
        print('Creating log file: %s' % log_file_name)
        _log_file = open(os.path.join(log_dir, log_file_name), 'w')
        
def close_log():
    
    global _log_file
    
    if log_file_open():
        _log_file.close()
        del _log_file
    
def restart_log():
    close_log()
    start_log()
    
def log_file_open():
    global _log_file
    return '_log_file' in globals()


### Date & time ###

def format_now():
    return datetime.now().strftime('%Y_%m_%d-%H:%M:%S')


### Iterators & collections ###

def to_chunks(iterable, chunk_size):
    
    chunk = []
    
    for element in iterable:
        
        chunk.append(element)
        
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
            
    if len(chunk) > 0:
        yield chunk


### argparse ###

def get_parser_bool_type(parser):

    def _bool_type(value):
        if isinstance(value, bool):
            return value
        if value.lower() in ['yes', 'true', 't', 'y', '1']:
            return True
        elif value.lower() in ['no', 'false', 'f', 'n', '0']:
            return False
        else:
            raise parser.error('"%s": unrecognized boolean value.' % value)
            
    return _bool_type

def get_parser_file_type(parser, must_exist = False):

    def _file_type(path):
    
        path = os.path.expanduser(path)
    
        if must_exist:
            if not os.path.exists(path):
                parser.error('File doesn\'t exist: %s' % path)
            elif not os.path.isfile(path):
                parser.error('Not a file: %s' % path)
            else:
                return path
        else:
        
            dir_path = os.path.dirname(path)
        
            if dir_path and not os.path.exists(dir_path):
                parser.error('Parent directory doesn\'t exist: %s' % dir_path)
            else:
                return path
    
    return _file_type

def get_parser_directory_type(parser, create_if_not_exists = False):
    
    def _directory_type(path):
    
        path = os.path.expanduser(path)
    
        if not os.path.exists(path):
            if create_if_not_exists:
            
                parent_path = os.path.dirname(path)
            
                if parent_path and not os.path.exists(parent_path):
                    parser.error('Cannot create empty directory (parent directory doesn\'t exist): %s' % path)
                else:
                    os.mkdir(path)
                    return path
            else:
                parser.error('Path doesn\'t exist: %s' % path)
        elif not os.path.isdir(path):
            parser.error('Not a directory: %s' % path)
        else:
            return path
        
    return _directory_type
