/*************************************************************************/
/*									 */
/*	Source code for use with Cubist Release 2.02			 */
/*	--------------------------------------------			 */
/*		   Copyright RuleQuest Research 2005			 */
/*									 */
/*	This code is provided "as is" without warranty of any kind,	 */
/*	either express or implied.  All use is at your own risk.	 */
/*									 */
/*************************************************************************/


/*************************************************************************/
/*									 */
/*		General data for Cubist					 */
/*									 */
/*************************************************************************/


	Attribute	ClassAtt=0,	/* attribute to use as class */
			LabelAtt=0;	/* attribute containing case label */

	int		MaxAtt,		/* max att number */
			MaxDiscrVal=3,	/* max discrete values for any att */
			Precision=2,	/* decimal places for target */
			MaxLabel=0,	/* max characters in case label */
			LineNo=0,	/* input line number */
			ErrMsgs=0,	/* errors found */
			AttExIn,	/* attribute exclusions/inclusions */
			TSBase=0,	/* base day for time stamps */
			Delimiter;	/* character at end of name */

	ItemNo		MaxItem=-1;	/* max data item number */

	Description	*Item;		/* data items */

	DiscrValue	*MaxAttVal,	/* number of values for each att */
			*Modal;		/* most frequent value for discr att */

	char		*SpecialStatus;	/* special att treatment */

	Definition	*AttDef;	/* definitions of implicit atts */

	String		Target,		/* name of dependent att */
		  	*AttName,	/* att names */
		  	**AttValName;	/* att value names */

	char		*IgnoredVals=0;	/* values of labels and ignored atts */
	int		IValsSize=0,	/* size of above */
			IValsOffset=0;	/* index of first free char */

	String		FileStem="undefined";
	char		Fn[512];	/* file name */

	FILE		*Mf=0;		/* file for reading models */

	float		*AttMean=Nil,	/* means of att values */
			*AttSD,		/* std dev ditto */
			Ceiling=1E38,	/* max allowable global prediction */
			Floor=-1E38;	/* min allowable global prediction */

	int		*AttPrec=Nil;	/* Attribute precision  */

	Description	*Instance=Nil,	/* training cases */
			Ref=Nil;	/* reference point */
	Index		KDTree=Nil;	/* index for searching training cases */
	ItemNo		MaxInstance=-1;	/* highest instance */
	float		*RSPredVal=Nil; /* tabulated RS predictions */
	NNEnvRec	GNNEnv;		/* global NN environment */

	unsigned char	*Tested;	/* used in BuildIndex */
	ItemCount	*ValFreq;	/* used in BuildIndex */

	Boolean		USEINSTANCES;
	float		EXTRAP=0.1,	/* allowed extrapolation from models */
			SAMPLE=0.0,	/* sample training proportion  */
			MAXD,		/* max distance for close neighbors */
			GlobalMean=0;
	int		MEMBERS=1,	/* models in committee */
			NN=5,		/* nearest neighbors to use */
			KRInit;
