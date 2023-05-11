#include "mlsa/MakeBlobPlacer.h"

#include <tcl.h>

#include "odb/db.h"
#include "mlsa/rtl_bp.h"
#include "ord/OpenRoad.hh"
#include "sta/StaMain.hh"

namespace sta {
    extern const char* mlsa_tcl_inits[];
}

extern "C" {
    extern int Mlsa_Init(Tcl_Interp* interp);
}

namespace ord {
    mlsa::BlobPlacer* makeBlobPlacer() {
        return new mlsa::BlobPlacer;
    }

    void initBlobPlacer(OpenRoad* openroad) {
        Tcl_Interp* tcl_interp = openroad->tclInterp();
        Mlsa_Init(tcl_interp);
        sta::evalTclInit(tcl_interp, sta::mlsa_tcl_inits);
        openroad->getBlobPlacer()->init(openroad->getDbNetwork(),
                                        openroad->getDb(),
                                        openroad->getSta(),
                                        openroad->getLogger(),
                                        openroad->getPartitionMgr());

    }

    void deleteBlobPlacer(mlsa::BlobPlacer* blob_placer) {
        delete blob_placer;
    }
} // namespace ord