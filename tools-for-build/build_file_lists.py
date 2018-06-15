#-*- mode: python; coding: utf-8-unix -*-

from cleavir_file_list import cleavir_file_list

SRC_CORE_FILES = \
    [
        'dummy',
        'clcenv',
        'mpPackage',
        'nativeVector',
        'environment',
        'activationFrame',
        'evaluator',
        'functor',
        'creator',
        'queue',
        'sharpEqualWrapper',
        'stacks',
        'weakKeyMapping',
        'weakHashTable',
        'weakPointer',
        'compiler',
        'genericFunction',
        'instance',
        'funcallableInstance',
        'cache',
        'float_to_string',
        'primitives',
        'random',
#        'cxxClass',
        'record',
        'debugger',
        'smallMap',
        'smallMultimap',
        'hashTable',
        'hashTableEq',
        'hashTableEql',
        'hashTableEqual',
        'hashTableEqualp',
        'numbers',
        'numerics',
        'num_arith',
        'numberToString',
        'num_co',
        'load',
        'bignum',
        'write_object',
        'write_array',
        'print',
        'sourceFileInfo',
        'symbolToEnumConverter',
        'core_globals',
        'externalObject',
        'myReadLine',
        'specialForm',
        'unixfsys',
        'lispList',
        'multiStringBuffer',
        'candoOpenMp',
        'foundation',
        'lambdaListHandler',
        'lispStream',
        'bits',
        'write_symbol',
        'corePackage',
        'lisp',
        'bundle',
        'write_ugly',
        'wrappedPointer',
        'serialize',
        'sexpLoadArchive',
        'sexpSaveArchive',
        'readtable',
        'float_to_digits',
        'pathname',
        'commandLineOptions',
        'exceptions',
        'commonLispUserPackage',
        'metaClass',
        'multipleValues',
        'testing',
        'predicates',
        'write_list',
        'package',
        'commonLispPackage',
        'allClSymbols',
        'keywordPackage',
        'extensionPackage',
        'array',
        'grayPackage',
        'closPackage',
        'cleavirPrimopsPackage',
        'cleavirEnvPackage',
        'compPackage',
        'bootStrapCoreSymbolMap',
        'cons',
        'symbol',
        'object',
        'arguments',
        'pointer',
        'iterator',
        'sysprop',
        'bformat',
        'backquote',
        'documentation',
        'lispReader',
        'singleDispatchGenericFunction',
        'singleDispatchMethod',
        'derivableCxxObject',
        'null',
        'character',
        'designators',
        'sequence',
        'loadTimeValues',
        'reader',
        'lightProfiler',
        'fileSystem',
        'intArray',
        'posixTime',
        'hwinfo',
        'clasp_ffi_package',
        'fli',
    ]

def collect_source_file(bld, path, name, extension = '.cc'):
    fullName = path + name
    if '.' not in name:
        fullName = fullName + extension
    node = bld.path.find_node(fullName)
    assert node != None, "Could not find %s" % (fullName)
    return node

def collect_c_source_files(bld, path, files, extension = '.cc'):
    result = []
    for name in files:
        result.append(collect_source_file(bld, path, name, extension = extension))
    return result

def collect_clasp_c_source_files(bld):
    result = collect_c_source_files(bld, 'src/gctools/', [
                 'gc_interface',
                 'boehmGarbageCollection',
                 'mpsGarbageCollection',
                 'hardErrors',
                 'source_info',
                 'threadlocal',
                 'gc_boot',
                 'interrupt',
                 'gcFunctions',
                 'gctoolsPackage',
                 'globals',
                 'gcStack',
                 'gcalloc',
                 'gcweak',
                 'memoryManagement',
                 'mygc.c']) + \
             collect_c_source_files(bld, 'src/clbind/', [
                 'adapter',
                 'class_rep',
                 'open',
                 'class_registry',
                 'link_compatibility',
                 'scope',
                 'inheritance',
                 'clbind',
                 'clbindPackage',
                 'class',
                 'derivable_class']) + \
             collect_c_source_files(bld, 'src/serveEvent/', [
                 'serveEvent',
                 'serveEventPackage']) + \
             collect_c_source_files(bld, 'src/sockets/', [
                 'sockets',
                 'socketsPackage']) + \
             collect_c_source_files(bld, 'src/llvmo/', [
                 'debugInfoExpose',
                 'debugLoc',
                 'llvmoDwarf',
                 'link_intrinsics',
                 'builtins',
                 'fastgf',
                 'intrinsics',
                 'insertPoint',
                 'irtests',
                 'llvmoExpose',
                 'llvmoPackage',
                 'clbindLlvmExpose']) + \
             collect_c_source_files(bld, 'src/asttooling/', [
                 'astVisitor',
                 'astExpose',
                 'clangTooling',
                 'asttoolingPackage',
                 'clangCompiler']) + \
             collect_c_source_files(bld, 'src/main/', ['main']) + \
             collect_c_source_files(bld, 'src/core/', SRC_CORE_FILES)

    return result
