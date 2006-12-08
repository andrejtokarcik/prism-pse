/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class simulator_SimulatorEngine */

#ifndef _Included_simulator_SimulatorEngine
#define _Included_simulator_SimulatorEngine
#ifdef __cplusplus
extern "C" {
#endif
#undef simulator_SimulatorEngine_ERROR
#define simulator_SimulatorEngine_ERROR -1L
#undef simulator_SimulatorEngine_OUTOFRANGE
#define simulator_SimulatorEngine_OUTOFRANGE -1L
#undef simulator_SimulatorEngine_NULL
#define simulator_SimulatorEngine_NULL 0L
#undef simulator_SimulatorEngine_NOT_LOADED
#define simulator_SimulatorEngine_NOT_LOADED 0L
#undef simulator_SimulatorEngine_PROBABILISTIC
#define simulator_SimulatorEngine_PROBABILISTIC 1L
#undef simulator_SimulatorEngine_NONDETERMINISTIC
#define simulator_SimulatorEngine_NONDETERMINISTIC 2L
#undef simulator_SimulatorEngine_STOCHASTIC
#define simulator_SimulatorEngine_STOCHASTIC 3L
#undef simulator_SimulatorEngine_UNDEFINED_INT
#define simulator_SimulatorEngine_UNDEFINED_INT -2147483647L
#undef simulator_SimulatorEngine_UNDEFINED_DOUBLE
#define simulator_SimulatorEngine_UNDEFINED_DOUBLE -1.0000000138484279E24
#undef simulator_SimulatorEngine_INFINITY
#define simulator_SimulatorEngine_INFINITY 1.0000000138484279E24
#undef simulator_SimulatorEngine_INTEGER
#define simulator_SimulatorEngine_INTEGER 1L
#undef simulator_SimulatorEngine_DOUBLE
#define simulator_SimulatorEngine_DOUBLE 2L
#undef simulator_SimulatorEngine_BOOLEAN
#define simulator_SimulatorEngine_BOOLEAN 3L
#undef simulator_SimulatorEngine_SIM_PATH_NUM_STEPS
#define simulator_SimulatorEngine_SIM_PATH_NUM_STEPS 0L
#undef simulator_SimulatorEngine_SIM_PATH_TIME
#define simulator_SimulatorEngine_SIM_PATH_TIME 1L
#undef simulator_SimulatorEngine_SIM_PATH_DEADLOCK
#define simulator_SimulatorEngine_SIM_PATH_DEADLOCK 2L
/*
 * Class:     simulator_SimulatorEngine
 * Method:    Set_Main_Log
 * Signature: (Lprism/PrismLog;)V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_Set_1Main_1Log
  (JNIEnv *, jclass, jobject);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    tidyUpEverything
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_tidyUpEverything
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    allocateStateSpace
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_allocateStateSpace
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    allocateModel
 * Signature: (III[I[III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_allocateModel
  (JNIEnv *, jclass, jint, jint, jint, jintArray, jintArray, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    setupAddTransition
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_setupAddTransition
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    setupAddStateReward
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_setupAddStateReward
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    setupAddTransitionReward
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_setupAddTransitionReward
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    allocatePath
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_allocatePath
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    defineVariable
 * Signature: (II)V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_defineVariable
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    startPath
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_startPath
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getPathSize
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getPathSize
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getPathData
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getPathData
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTimeSpentInPathState
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTimeSpentInPathState
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getCumulativeTimeSpentInPathState
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getCumulativeTimeSpentInPathState
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getStateRewardOfPathState
 * Signature: (II)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getStateRewardOfPathState
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTransitionRewardOfPathState
 * Signature: (II)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTransitionRewardOfPathState
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalStateRewardOfPathState
 * Signature: (II)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalStateRewardOfPathState
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalTransitionRewardOfPathState
 * Signature: (II)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalTransitionRewardOfPathState
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalPathTime
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalPathTime
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalPathReward
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalPathReward
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalStateReward
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalStateReward
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getTotalTransitionReward
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getTotalTransitionReward
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    isPathLooping
 * Signature: ()Z
 */
JNIEXPORT jboolean JNICALL Java_simulator_SimulatorEngine_isPathLooping
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loopStart
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loopStart
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loopEnd
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loopEnd
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getChosenIndexOfOldUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getChosenIndexOfOldUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    makeManualUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_makeManualUpdate__I
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    makeManualUpdate
 * Signature: (ID)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_makeManualUpdate__ID
  (JNIEnv *, jclass, jint, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    doAutomaticChoices
 * Signature: (IZ)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_doAutomaticChoices
  (JNIEnv *, jclass, jint, jboolean);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    doBacktrack
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_doBacktrack
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    doRemovePrecedingStates
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_doRemovePrecedingStates
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    calculateOldUpdates
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_calculateOldUpdates
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    finishedWithOldUpdates
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_finishedWithOldUpdates
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getNumUpdates
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getNumUpdates
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getActionIndexOfUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getActionIndexOfUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getModuleIndexOfUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getModuleIndexOfUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getProbabilityOfUpdate
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getProbabilityOfUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getNumAssignmentsOfUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getNumAssignmentsOfUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getAssignmentVariableIndexOfUpdate
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getAssignmentVariableIndexOfUpdate
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getAssignmentValueOfUpdate
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getAssignmentValueOfUpdate
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getDistributionIndexOfUpdate
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getDistributionIndexOfUpdate
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    allocatePCTLManager
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_allocatePCTLManager
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    allocateSampling
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_allocateSampling
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    exportBinary
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_exportBinary
  (JNIEnv *, jclass, jstring);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    doSampling
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_doSampling
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getSamplingResult
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_getSamplingResult
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getNumReachedMaxPath
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_getNumReachedMaxPath
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    stopSampling
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_stopSampling
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadProposition
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadProposition
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryProposition
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryProposition__I
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryProposition
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryProposition__II
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryIsInitial
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryIsInitial__
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryIsInitial
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryIsInitial__I
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryIsDeadlock
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryIsDeadlock__
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryIsDeadlock
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryIsDeadlock__I
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    findPathFormulaIndex
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_findPathFormulaIndex
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryPathFormula
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_queryPathFormula
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    queryPathFormulaNumeric
 * Signature: (I)D
 */
JNIEXPORT jdouble JNICALL Java_simulator_SimulatorEngine_queryPathFormulaNumeric
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalConstant
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalConstant
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealConstant
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealConstant
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createIntegerVar
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createIntegerVar
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createBooleanVar
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createBooleanVar
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createDouble
 * Signature: (D)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createDouble
  (JNIEnv *, jclass, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createInteger
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createInteger
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createBoolean
 * Signature: (Z)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createBoolean
  (JNIEnv *, jclass, jboolean);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createCeil
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createCeil
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createFloor
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createFloor
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalPow
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalPow
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealPow
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealPow
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createMod
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createMod
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNot
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNot
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createAnd
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createAnd
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createOr
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createOr
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalMax
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalMax
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalMin
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalMin
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealMax
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealMax
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealMin
 * Signature: ([I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealMin
  (JNIEnv *, jclass, jintArray);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalTimes
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalTimes
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalPlus
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalPlus
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalMinus
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalMinus
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealTimes
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealTimes
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealPlus
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealPlus
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealMinus
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealMinus
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createDivide
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createDivide
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createIte
 * Signature: (III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createIte
  (JNIEnv *, jclass, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealIte
 * Signature: (III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealIte
  (JNIEnv *, jclass, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalEquals
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalEquals
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealEquals
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealEquals
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalNotEquals
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalNotEquals
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealNotEquals
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealNotEquals
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalLessThan
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalLessThan
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealLessThan
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealLessThan
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalGreaterThan
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalGreaterThan
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealGreaterThan
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealGreaterThan
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalLessThanEqual
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalLessThanEqual
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealLessThanEqual
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealLessThanEqual
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createNormalGreaterThanEqual
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createNormalGreaterThanEqual
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createRealGreaterThanEqual
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createRealGreaterThanEqual
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlBoundedUntil
 * Signature: (IIDD)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlBoundedUntil
  (JNIEnv *, jclass, jint, jint, jdouble, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlBoundedUntilNegated
 * Signature: (IIDD)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlBoundedUntilNegated
  (JNIEnv *, jclass, jint, jint, jdouble, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlUntil
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlUntil
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlUntilNegated
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlUntilNegated
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlNext
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlNext
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlReachability
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlReachability
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlCumulative
 * Signature: (ID)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlCumulative
  (JNIEnv *, jclass, jint, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadPctlInstantanious
 * Signature: (ID)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadPctlInstantanious
  (JNIEnv *, jclass, jint, jdouble);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadProbQuestion
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadProbQuestion
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    loadRewardQuestion
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_loadRewardQuestion
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createCommand
 * Signature: (IIII)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createCommand
  (JNIEnv *, jclass, jint, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    addUpdate
 * Signature: (III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_addUpdate
  (JNIEnv *, jclass, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    addAssignment
 * Signature: (III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_addAssignment
  (JNIEnv *, jclass, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createStateReward
 * Signature: (II)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createStateReward
  (JNIEnv *, jclass, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    createTransitionReward
 * Signature: (III)I
 */
JNIEXPORT jint JNICALL Java_simulator_SimulatorEngine_createTransitionReward
  (JNIEnv *, jclass, jint, jint, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    printExpression
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_printExpression
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    expressionToString
 * Signature: (I)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_simulator_SimulatorEngine_expressionToString
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    deleteExpression
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_deleteExpression
  (JNIEnv *, jclass, jint);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    modelToString
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_simulator_SimulatorEngine_modelToString
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    pathToString
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_simulator_SimulatorEngine_pathToString
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    printCurrentUpdates
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_simulator_SimulatorEngine_printCurrentUpdates
  (JNIEnv *, jclass);

/*
 * Class:     simulator_SimulatorEngine
 * Method:    getLastErrorMessage
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_simulator_SimulatorEngine_getLastErrorMessage
  (JNIEnv *, jclass);

#ifdef __cplusplus
}
#endif
#endif
