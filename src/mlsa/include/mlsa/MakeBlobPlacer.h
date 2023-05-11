#pragma once

namespace mlsa {
    class BlobPlacer;
}

namespace ord {
    class OpenRoad;

    mlsa::BlobPlacer* makeBlobPlacer();

    void initBlobPlacer(OpenRoad* openroad);

    void deleteBlobPlacer(mlsa::BlobPlacer* blob_placer);
} // namespace ord