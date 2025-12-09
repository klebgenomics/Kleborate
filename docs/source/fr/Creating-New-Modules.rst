################################################
Création de nouveaux modules
################################################


#. 
   **Créer un nouveau répertoire de modules**\ :


   * Naviguez dans le répertoire ``modules`` du projet Kleborate.
   * Créez un nouveau répertoire pour votre module. Le nom du répertoire doit décrire la fonctionnalité du module.

#. 
   **Créer le fichier Python du module**\ :


   * Dans le répertoire nouvellement créé, créez un fichier Python avec le même nom que le répertoire. Ce fichier contiendra l'implémentation du module.
   * Par exemple, si votre module est nommé ``enterobacterales_species``\ , créez un fichier nommé ``enterobacterales_species.py`` dans le répertoire ``enterobacterales_species``.

#. 
   **Implémenter la fonctionnalité du module**\ :


   * Définissez les fonctionnalités de votre module dans le fichier Python.
   * Vous devez définir les fonctions pour ajouter des options CLI (command line), vérifier les options CLI, vérifier les dépendances externes, obtenir les en-têtes de module, obtenir les résultats de module, et toute autre fonctionnalité que votre module exige.
   * Assurez-vous que votre module respecte les conventions et lignes directrices des modules Kleborate existants.

#. 
   **Ajouter les options CLI**\ :


   * Implémentez une fonction pour ajouter des options de ligne de commande spécifiques à votre module. Cette fonction devrait accepter un objet ``argparse.ArgumentParser`` comme argument et ajouter des options en utilisant ses méthodes.
   * La fonction doit être nommée ``add_cli_options(parser)`` et retourner le groupe d'arguments créé pour les options du module.

#. 
   **Check CLI Options**\ :


   * Implémenter une fonction pour vérifier les options CLI fournies par l'utilisateur. S'assurer que les options fournies sont valides et répondent aux exigences du module.
   * Cette fonction doit être nommée ``check_cli_options(args)`` et soulever les erreurs ou les avertissements appropriés si les options sont invalides.

#. 
   **Vérifier les dépendances du programme externe**\ :


   * Implémentez une fonction pour vérifier les dépendances externes requises par votre module.
   * Cette fonction doit être nommée ``check_external_programs()`` et retourner un ensemble contenant les noms de programmes externes requis par le module.

#. 
   **Définir les modules prérequis**\ :


   * Définissez une fonction nommée ``prerequisite_modules()`` qui renvoie une liste de noms de modules dont votre module dépend.
   * Cela garantit que les modules préalables sont exécutés avant votre module.

#. 
   **Définir les en-têtes de module**\ :


   * Implémenter une fonction pour définir les en-têtes de sortie du module.
   * Cette fonction devrait renvoyer deux listes : une pour les en-têtes complets et une pour les en-têtes stdout.

#. 
   **Définir les résultats du module**\ :


   * Implémentez une fonction pour obtenir les résultats produits par votre module.
   * Cette fonction doit accepter les arguments nécessaires comme l'assemblage, l'index minimap2, les arguments de ligne de commande et les autres données requises.
   * Il doit renvoyer un dictionnaire contenant les résultats.

#. 
   **Testez votre module**\ :


   * Testez votre module avec différentes entrées et différents scénarios pour s'assurer qu'il se comporte comme prévu.
   * Vérifier que le module produit des résultats corrects et gère les erreurs.

#. 
   **Documenter votre module**\ :


   * Écrivez la documentation pour votre module, y compris son but, son utilisation, ses options et toute autre information pertinente.
   * Mettre à jour la documentation principale de Kleborate pour inclure des informations sur votre module.
