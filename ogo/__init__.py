'''Like ogopogo, osteoporosis is elusive'''

__all__ = (
  '__version__',
  'version_info',
)
try:
  from pbr.version import VersionInfo
except ImportError:
  __version__ = '0+unknown'
  version_info = (0,)
else:
  _v = VersionInfo('ogo').semantic_version()
  __version__ = _v.release_string()
  version_info = _v.version_tuple()
