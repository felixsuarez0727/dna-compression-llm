from dna_compression.detection.local_models import _place_model_on_available_device


class FakeCuda:
    def __init__(self, available):
        self.available = available

    def is_available(self):
        return self.available


class FakeTorch:
    def __init__(self, cuda_available):
        self.cuda = FakeCuda(cuda_available)


class FakeModel:
    def __init__(self, cuda_error=None):
        self.cuda_error = cuda_error
        self.devices = []

    def to(self, device):
        self.devices.append(device)
        if device == "cuda" and self.cuda_error:
            raise self.cuda_error
        return self


def test_places_model_on_cuda_when_available():
    logs = []
    model = FakeModel()

    placed_model, device = _place_model_on_available_device(
        model, FakeTorch(cuda_available=True), logs.append, "Test model"
    )

    assert placed_model is model
    assert device == "cuda"
    assert model.devices == ["cuda"]
    assert logs == ["  Test model loaded on CUDA."]


def test_places_model_on_cpu_when_cuda_is_unavailable():
    logs = []
    model = FakeModel()

    placed_model, device = _place_model_on_available_device(
        model, FakeTorch(cuda_available=False), logs.append, "Test model"
    )

    assert placed_model is model
    assert device == "cpu"
    assert model.devices == ["cpu"]
    assert logs == ["  Test model loaded on CPU."]


def test_falls_back_to_cpu_when_cuda_placement_fails():
    logs = []
    model = FakeModel(cuda_error=RuntimeError("out of memory"))

    placed_model, device = _place_model_on_available_device(
        model, FakeTorch(cuda_available=True), logs.append, "Test model"
    )

    assert placed_model is model
    assert device == "cpu"
    assert model.devices == ["cuda", "cpu"]
    assert logs == [
        "  Test model: CUDA placement failed (out of memory); falling back to CPU.",
        "  Test model loaded on CPU.",
    ]