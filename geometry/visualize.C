#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>

void visualize() {
    // 1. Inicializar el TGeoManager primero para evitar el error de puntero nulo
    if (!gGeoManager) {
        new TGeoManager("geom", "ATTPC Geometry Visualizer");
    }

    // 2. Abrir el archivo ROOT existente
    TFile *f = new TFile("ATTPC_H300torr_RT.root", "READ");
    if (f->IsZombie()) {
        std::cerr << "Error: No se pudo abrir el archivo ATTPC_H300torr_RT.root" << std::endl;
        return;
    }

    // 3. Extraer el volumen TOP
    TGeoVolume *TOP = nullptr;
    f->GetObject("TOP", TOP);

    if (TOP) {
        // 4. Configurar el gestor de geometría y cerrar
        gGeoManager->SetTopVolume(TOP);
        gGeoManager->CloseGeometry();

        // 5. Mostrar la visualización en 3D (OpenGL)
        TOP->Draw("ogl");
        std::cout << "¡Geometría cargada y lista en el visor 3D!" << std::endl;
    } else {
        std::cerr << "Error: No se encontró el volumen 'TOP' dentro del archivo." << std::endl;
    }
}
