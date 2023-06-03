final int DIMENSION = 600;                                                               // Taille de la grille
final double convolution[][] = {{0.05, 0.2, 0.05}, {0.2, -1., 0.2}, {0.05, 0.2, 0.05}}; // Matrice de convolution pour le Laplacien

double[][] concentrationA = new double[DIMENSION][DIMENSION];                          // Concentration en A
double[][] concentrationB = new double[DIMENSION][DIMENSION];                          // Concentration en B
double[][] variationA = new double[DIMENSION][DIMENSION];                              // Variation de la concentration en A
double[][] variationB = new double[DIMENSION][DIMENSION];                              // Variation de la concentration en B

double diffusionA = 1.0;
double diffusionB = 0.5;
double feedRate = 0.055;
double killRate = 0.062;

void setup() {
  size(600, 600);
  smooth();
  colorMode(RGB);
  initializeState();
}

void draw() {
  for (int p = 0; p < 100; p++) {
    reactionDiffusion();
  }

  background(255);                                                               // Fond blanc

  for (int i = 0; i < DIMENSION; i++) {
    for (int j = 0; j < DIMENSION; j++) {
      float colorValue = (float) (255 - concentrationA[i][j] * 255);              // Calcul de la valeur de gris en inversant la concentration
      stroke(colorValue);                                                        // Couleur de trait en gris
      point(i, j);
    }
  }
}

void initializeState() {
  for (int i = 0; i < DIMENSION; i++) {
    for (int j = 0; j < DIMENSION; j++) {
      concentrationA[i][j] = 1.0;
      concentrationB[i][j] = 0.0;
      variationA[i][j] = variationB[i][j] = 0;
    }
  }

  // Initialisation d'une petite zone ensemencée avec B=1
  int seedSize = 30;
  int seedX = DIMENSION / 2 - seedSize / 2;
  int seedY = DIMENSION / 2 - seedSize / 2;
  for (int i = 0; i < seedSize; i++) {
    for (int j = 0; j < seedSize; j++) {
      concentrationB[seedX + i][seedY + j] = 1.0;
    }
  }
}

void reactionDiffusion() {
  for (int i = 0; i < DIMENSION; i++) {
    for (int j = 0; j < DIMENSION; j++) {
      double a = concentrationA[i][j];
      double b = concentrationB[i][j];
      double reactionA = diffusionA * laplacian(i, j, concentrationA) - a * b * b + feedRate * (1 - a);
      double reactionB = diffusionB * laplacian(i, j, concentrationB) + a * b * b - (killRate + feedRate) * b;
      variationA[i][j] = reactionA;
      variationB[i][j] = reactionB;
    }
  }

  for (int i = 0; i < DIMENSION; i++) {
    for (int j = 0; j < DIMENSION; j++) {
      concentrationA[i][j] += variationA[i][j];
      concentrationB[i][j] += variationB[i][j];
      concentrationA[i][j] = constrain((float) concentrationA[i][j], 0, 1);
      concentrationB[i][j] = constrain((float) concentrationB[i][j], 0, 1);
    }
  }
}

double laplacian(int x, int y, double[][] grid) {
  double laplace = 0;
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      int xi = (x + i + DIMENSION) % DIMENSION;                                    // Gestion des bords périodiques
      int yj = (y + j + DIMENSION) % DIMENSION;                                    // Gestion des bords périodiques
      laplace += convolution[i + 1][j + 1] * grid[xi][yj];                         // Calcul de la somme pondérée des concentrations environnantes
    }
  }
  return laplace;
}
