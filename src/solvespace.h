//-----------------------------------------------------------------------------
// All declarations not grouped specially elsewhere.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------

#ifndef SOLVESPACE_H
#define SOLVESPACE_H

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <memory>

#define EIGEN_NO_DEBUG
#undef Success
#include <Eigen/SparseCore>

#include "dsc.h"
#include "polygon.h"
#include "srf/surface.h"
#include "render/render.h"
#include "expr.h"
#include "sketch.h"
#include "ttf.h"
#include "ui.h"

#include "platform/platform.h"

namespace SolveSpace {

using std::min;
using std::max;
using std::swap;
using std::fabs;

enum class Unit : uint32_t {
    MM = 0,
    INCHES,
    METERS,
    FEET_INCHES
};

class Pixmap;

enum class SolveResult : uint32_t {
    OKAY                     = 0,
    DIDNT_CONVERGE           = 10,
    REDUNDANT_OKAY           = 11,
    REDUNDANT_DIDNT_CONVERGE = 12,
    TOO_MANY_UNKNOWNS        = 20
};

// Utility functions that are provided in the platform-independent code.
class utf8_iterator {
    const char *p, *n;
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = char32_t;
    using difference_type = std::ptrdiff_t;
    using pointer = char32_t*;
    using reference = char32_t&;

    utf8_iterator(const char *p) : p(p), n(NULL) {}
    bool           operator==(const utf8_iterator &i) const { return p==i.p; }
    bool           operator!=(const utf8_iterator &i) const { return p!=i.p; }
    ptrdiff_t      operator- (const utf8_iterator &i) const { return p -i.p; }
    utf8_iterator& operator++()    { **this; p=n; n=NULL; return *this; }
    utf8_iterator  operator++(int) { utf8_iterator t(*this); operator++(); return t; }
    char32_t       operator*();
    const char*    ptr() const { return p; }
};
class ReadUTF8 {
    const std::string &str;
public:
    ReadUTF8(const std::string &str) : str(str) {}
    utf8_iterator begin() const { return utf8_iterator(&str[0]); }
    utf8_iterator end()   const { return utf8_iterator(&str[0] + str.length()); }
};









class SolveSpaceUI {
public:
    TextWindow                 *pTW;
    TextWindow                 &TW;
    GraphicsWindow              GW;

    // The state for undo/redo
    typedef struct UndoState {
        IdList<Group,hGroup>            group;
        List<hGroup>                    groupOrder;
        IdList<Request,hRequest>        request;
        IdList<Constraint,hConstraint>  constraint;
        ParamList                       param;
        IdList<Style,hStyle>            style;
        hGroup                          activeGroup;

        void Clear() {
            group.Clear();
            request.Clear();
            constraint.Clear();
            param.Clear();
            style.Clear();
        }
    } UndoState;
    enum { MAX_UNDO = 100 };
    typedef struct {
        UndoState   d[MAX_UNDO];
        int         cnt;
        int         write;
    } UndoStack;
    UndoStack   undo;
    UndoStack   redo;

    std::map<Platform::Path, std::shared_ptr<Pixmap>, Platform::PathLess> images;
    bool ReloadLinkedImage(const Platform::Path &saveFile, Platform::Path *filename,
                           bool canCancel);

    void UndoEnableMenus();
    void UndoRemember();
    void UndoUndo();
    void UndoRedo();
    void PushFromCurrentOnto(UndoStack *uk);
    void PopOntoCurrentFrom(UndoStack *uk);
    void UndoClearState(UndoState *ut);
    void UndoClearStack(UndoStack *uk);

    // Little bits of extra configuration state
    enum { MODEL_COLORS = 8 };
    RgbaColor modelColor[MODEL_COLORS];
    Vector   lightDir[2];
    double   lightIntensity[2];
    double   ambientIntensity;
    double   chordTol;
    double   chordTolCalculated;
    int      maxSegments;
    double   exportChordTol;
    int      exportMaxSegments;
    int      timeoutRedundantConstr; //milliseconds
    int      animationSpeed; //milliseconds
    double   cameraTangent;
    double   gridSpacing;
    double   exportScale;
    double   exportOffset;
    bool     arcDimDefaultDiameter;
    bool     showFullFilePath;
    bool     fixExportColors;
    bool     exportBackgroundColor;
    bool     drawBackFaces;
    bool     showContourAreas;
    bool     checkClosedContour;
    bool     cameraNav;
    bool     turntableNav;
    bool     immediatelyEditDimension;
    bool     automaticLineConstraints;
    bool     showToolbar;
    Platform::Path screenshotFile;
    RgbaColor backgroundColor;
    bool     exportShadedTriangles;
    bool     exportPwlCurves;
    bool     exportCanvasSizeAuto;
    bool     exportMode;
    struct {
        double  left;
        double  right;
        double  bottom;
        double  top;
    }        exportMargin;
    struct {
        double  width;
        double  height;
        double  dx;
        double  dy;
    }        exportCanvas;
    struct {
        double  depth;
        double  safeHeight;
        int     passes;
        double  feed;
        double  plungeFeed;
    }        gCode;

    Unit     viewUnits;
    int      afterDecimalMm;
    int      afterDecimalInch;
    int      afterDecimalDegree;
    bool     useSIPrefixes;
    int      autosaveInterval; // in minutes
    bool     explode;
    double   explodeDistance;

    std::string MmToString(double v, bool editable=false);
    std::string MmToStringSI(double v, int dim = 0);
    std::string DegreeToString(double v);
    double ExprToMm(Expr *e);
    double StringToMm(const std::string &s);
    const char *UnitName();
    double MmPerUnit();
    int UnitDigitsAfterDecimal();
    void SetUnitDigitsAfterDecimal(int v);
    double ChordTolMm();
    double ExportChordTolMm();
    int GetMaxSegments();
    bool usePerspectiveProj;
    double CameraTangent();

    // Some stuff relating to the tangent arcs created non-parametrically
    // as special requests.
    double tangentArcRadius;
    bool tangentArcManual;
    bool tangentArcModify;

    // The platform-dependent code calls this before entering the msg loop
    void Init();
    void Exit();

    // File load/save routines, including the additional files that get
    // loaded when we have link groups.
    FILE        *fh;
    void AfterNewFile();
    void AddToRecentList(const Platform::Path &filename);
    Platform::Path saveFile;
    bool        fileLoadError;
    bool        unsaved;
    typedef struct {
        char        type;
        const char *desc;
        char        fmt;
        void       *ptr;
    } SaveTable;
    static const SaveTable SAVED[];
    void SaveUsingTable(const Platform::Path &filename, int type);
    void LoadUsingTable(const Platform::Path &filename, char *key, char *val);
    struct {
        Group        g;
        Request      r;
        Entity       e;
        Param        p;
        Constraint   c;
        Style        s;
    } sv;
    static void MenuFile(Command id);
    void Autosave();
    void RemoveAutosave();
    static constexpr size_t MAX_RECENT = 8;
    static constexpr const char *SKETCH_EXT = "slvs";
    static constexpr const char *BACKUP_EXT = "slvs~";
    std::vector<Platform::Path> recentFiles;
    bool Load(const Platform::Path &filename);
    bool GetFilenameAndSave(bool saveAs);
    bool OkayToStartNewFile();
    hGroup CreateDefaultDrawingGroup();
    void UpdateWindowTitles();
    void ClearExisting();
    void NewFile();
    bool SaveToFile(const Platform::Path &filename);
    bool LoadAutosaveFor(const Platform::Path &filename);
    std::function<void(const Platform::Path &filename, bool is_saveAs, bool is_autosave)> OnSaveFinished;
    bool LoadFromFile(const Platform::Path &filename, bool canCancel = false);
    void UpgradeLegacyData();
    bool LoadEntitiesFromFile(const Platform::Path &filename, EntityList *le,
                              SMesh *m, SShell *sh);
    bool LoadEntitiesFromSlvs(const Platform::Path &filename, EntityList *le,
                              SMesh *m, SShell *sh);
    bool ReloadAllLinked(const Platform::Path &filename, bool canCancel = false);
    // And the various export options
    void ExportAsPngTo(const Platform::Path &filename);
    void ExportMeshTo(const Platform::Path &filename);
    void ExportMeshAsStlTo(FILE *f, SMesh *sm);
    void ExportMeshAsObjTo(FILE *fObj, FILE *fMtl, SMesh *sm);
    void ExportMeshAsThreeJsTo(FILE *f, const Platform::Path &filename,
                               SMesh *sm, SOutlineList *sol);
    void ExportMeshAsVrmlTo(FILE *f, const Platform::Path &filename, SMesh *sm);
    void ExportViewOrWireframeTo(const Platform::Path &filename, bool exportWireframe);
    void ExportSectionTo(const Platform::Path &filename);
    void ExportWireframeCurves(SEdgeList *sel, SBezierList *sbl,
                               VectorFileWriter *out);
    void ExportLinesAndMesh(SEdgeList *sel, SBezierList *sbl, SMesh *sm,
                            Vector u, Vector v,
                            Vector n, Vector origin,
                            double cameraTan,
                            VectorFileWriter *out);

    static void MenuAnalyze(Command id);

    // Additional display stuff
    struct {
        SContour    path;
        hEntity     point;
    } traced;
    SEdgeList nakedEdges;
    struct {
        bool        draw;
        Vector      ptA;
        Vector      ptB;
    } extraLine;
    struct {
        bool        draw, showOrigin;
        Vector      pt, u, v;
    } justExportedInfo;
    struct {
        bool   draw;
        bool   dirty;
        Vector position;
    } centerOfMass;

    class Clipboard {
    public:
        List<ClipboardRequest>  r;
        List<Constraint>        c;

        void Clear();
        bool ContainsEntity(hEntity old);
        hEntity NewEntityFor(hEntity old);
    };
    Clipboard clipboard;

    void MarkGroupDirty(hGroup hg, bool onlyThis = false);
    void MarkGroupDirtyByEntity(hEntity he);

    // Consistency checking on the sketch: stuff with missing dependencies
    // will get deleted automatically.
    struct {
        int     requests;
        int     groups;
        int     constraints;
        int     nonTrivialConstraints;
    } deleted;
    bool GroupExists(hGroup hg);
    bool EntityExists(hEntity he);
    bool GroupsInOrder(hGroup hbefore, hGroup hafter);
    bool PruneGroup(Group *g);
    static void ShowNakedEdges(bool reportOnlyWhenNotOkay);

    enum class Generate : uint32_t {
        DIRTY,
        ALL,
        REGEN,
        UNTIL_ACTIVE,
    };

    void GenerateAll(Generate type = Generate::DIRTY, bool andFindFree = false);
    void SolveGroup(hGroup hg, bool andFindFree);
    void SolveGroupAndReport(hGroup hg, bool andFindFree);
    SolveResult TestRankForGroup(hGroup hg, int *rank = NULL);
    void WriteEqSystemForGroup(hGroup hg);
    void MarkDraggedParams();
    void ForceReferences();
    void UpdateCenterOfMass();

    bool ActiveGroupsOkay();

    // The system to be solved.
    System     *pSys;
    System     &sys;

    // All the TrueType fonts in memory
    TtfFontList fonts;

    // Everything has been pruned, so we know there's no dangling references
    // to entities that don't exist. Before that, we mustn't try to display
    // the sketch!
    bool allConsistent;

    bool scheduledGenerateAll;
    bool scheduledShowTW;
    Platform::TimerRef refreshTimer;
    Platform::TimerRef autosaveTimer;
    void Refresh();
    void ScheduleShowTW();
    void ScheduleGenerateAll();
    void ScheduleAutosave();

    static void MenuHelp(Command id);

    void Clear();

    // We allocate TW and sys on the heap to work around an MSVC problem
    // where it puts zero-initialized global data in the binary (~30M of zeroes)
    // in release builds.
    SolveSpaceUI()
        : pTW(new TextWindow()), TW(*pTW),
          pSys(new System()), sys(*pSys) {}

    ~SolveSpaceUI() {
        delete pTW;
        delete pSys;
    }
};

void ImportDxf(const Platform::Path &file);
void ImportDwg(const Platform::Path &file);
bool LinkIDF(const Platform::Path &filename, EntityList *le, SMesh *m, SShell *sh);
bool LinkStl(const Platform::Path &filename, EntityList *le, SMesh *m, SShell *sh);

extern SolveSpaceUI SS;
extern Sketch SK;

} // namespace SolveSpace

#endif
