"""Public failure types for NTv3 planning and artifact generation."""


class Ntv3Error(Exception):
    """Base class for expected, user-actionable NTv3 failures."""


class RequestError(Ntv3Error, ValueError):
    """Raised when a request is incomplete, inconsistent, or ambiguous."""


class ArtifactError(Ntv3Error, ValueError):
    """Raised when required run artifacts are missing or inconsistent."""


class ExternalServiceError(Ntv3Error, RuntimeError):
    """Raised when a public data service cannot satisfy a request."""
